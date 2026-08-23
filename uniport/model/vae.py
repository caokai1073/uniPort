'''
Author: Kai Cao
Modified from SCALEX
'''

import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from tqdm import tqdm
from collections import defaultdict
from .layer import Encoder, Decoder
from .loss import kl_div, unbalanced_ot


class VAE(nn.Module):
    """
    Variational Autoencoder framework
    """
    def __init__(self, enc, dec, ref_id, n_domain, mode):
        """
        Parameters
        ----------
        enc
            Encoder structure config
        dec
            Decoder structure config
        ref_id
            ID of reference dataset
        n_domain
            The number of different domains
        mode
            Choose from ['h', 'v', 'd']

        """
        super().__init__()


        x_dim = {}

        for key in dec.keys():
            x_dim[key] = dec[key][-1][1]
        self.z_dim = enc[-1][1]

        self.encoder = Encoder(x_dim, enc, mode)
        self.decoder = Decoder(self.z_dim, dec)

        self.n_domain = n_domain
        self.ref_id = ref_id
  
    def load_model(self, path):
        """
        Load trained model parameters dictionary.
        Parameters
        ----------
        path
            file path that stores the model parameters
        """
        pretrained_dict = torch.load(path, map_location='cpu', weights_only=True)
        model_dict = self.state_dict()
        pretrained_dict = {k: v for k, v in pretrained_dict.items() if k in model_dict}
        model_dict.update(pretrained_dict) 
        self.load_state_dict(model_dict)
        
    def encodeBatch(
            self, 
            dataloader, 
            num_gene,
            mode='h',  
            out='latent', 
            batch_id=0,
            pred_id=1,
            device='cuda', 
            eval=False,
        ):
        """
        Inference

        Parameters
        ----------
        dataloader
            An iterable over the given dataset for inference.
        num_gene
            List of number of genes in different datasets. 
        mode
            Choose from ['h', 'v', 'd']
            If 'h', integrate data with common genes
            If 'v', integrate data profiled from the same cells
            If 'd', inetrgate data without common genes
            Default: 'h'.
        out
            Output of uniPort. Choose from ['latent', 'project', 'predict']. 
            If out='latent', train the network and output cell embeddings. 
            If out='project', project data into the latent space and output cell embeddings. 
            If out='predict', project data into the latent space and output cell embeddings through a specified decoder. 
            Default: 'latent'.
        batch_id
            Choose which encoder to project data when mode='d'. Default: 0.
        pred_id
            Choose which decoder to reconstruct data when out='predict'.
        device
            'cuda' or 'cpu' for . Default: 'cuda'.
        eval
            If True, set the model to evaluation mode. If False, set the model to train mode. Default: False.
        
        Returns
        -------
        output
            Cell embeddings (if out='latent' or 'project') or Predicted data (if out='predict').
        """

        self.to(device)
        if eval:
            self.eval()
        else:
            self.train()

        # Inference never backpropagates: the results were already detached
        # before being collected, so building the graph was pure overhead.
        with torch.no_grad():
            output = []
            if out == 'latent' or out == 'project':

                if mode == 'v':
                    for x, idx in dataloader:
                        x = x.to(device, non_blocking=True).float()
                        z = self.encoder(x[:, 0:num_gene[0]].contiguous(), 0)[1]
                        output.append(z.cpu())
                    output = torch.cat(output).numpy()

                elif mode == 'd':
                    for x, y, idx in dataloader:
                        x = x[:, 0:num_gene[batch_id]].contiguous().to(device, non_blocking=True).float()
                        y = y.to(device, non_blocking=True).long()
                        loc = torch.where(y==batch_id)[0]
                        x = x[loc]
                        z = self.encoder(x, batch_id)[1] # z, mu, var
                        output.append(z.cpu())
                    output = torch.cat(output).numpy()

                elif mode == 'h':
                    output = np.zeros((dataloader.dataset.shape[0], self.z_dim))
                    for x, y, idx in dataloader:
                        x_c = x[:, 0:num_gene[self.n_domain]].contiguous().to(device, non_blocking=True).float()
                        z = self.encoder(x_c, 0)[1]
                        output[idx.numpy()] = z.cpu().numpy()

            elif out == 'predict':

                if mode == 'v':
                    for x, idx in dataloader:
                        x = x.to(device, non_blocking=True).float()
                        z = self.encoder(x[:, 0:num_gene[0]].contiguous(), 0)[1]
                        recon = self.decoder(z, pred_id)
                        output.append(recon.cpu())
                    output = torch.cat(output).numpy()

                elif mode == 'd':
                    for x, y, idx in dataloader:
                        x = x[:, 0:num_gene[batch_id]].contiguous().to(device, non_blocking=True).float()
                        y = y.to(device, non_blocking=True).long()
                        loc = torch.where(y==batch_id)[0]
                        x = x[loc]
                        z = self.encoder(x, batch_id)[1] # z, mu, var
                        recon = self.decoder(z, pred_id)
                        output.append(recon.cpu())
                    output = torch.cat(output).numpy()

                elif mode == 'h':
                    output = np.zeros((dataloader.dataset.shape[0], num_gene[pred_id]))
                    for x, y, idx in dataloader:
                        x_c = x[:, 0:num_gene[self.n_domain]].contiguous().to(device, non_blocking=True).float()
                        z = self.encoder(x_c, 0)[1]
                        recon = self.decoder(z, pred_id+1)
                        output[idx.numpy()] = recon.cpu().numpy()

        return output

    def fit(
            self, 
            dataloader, 
            tran,
            num_cell,
            num_gene,
            mode='h',
            loss_type='BCE',
            label_weight=None,
            Prior=None,
            save_OT=False,
            use_specific=True,
            lambda_s=0.5,
            lambda_kl=0.5,
            lambda_recon=1.0,
            lambda_ot=1.0,
            reg=0.1,
            reg_m=1.0,
            lr=2e-4,
            max_iteration=30000,
            early_stopping=None,
            device='cuda:0',  
            verbose=False,
        ):
        """
        train VAE

        Parameters
        ----------
        dataloader
            An iterable over the given dataset for training.
        tran
            A global OT plan. tran={} if save_OT=False in function.py.
        num_cell
            List of number of cells in different datasets. 
        num_gene
            List of number of genes in different datasets. 
        mode
            Choose from ['h', 'v', 'd']
            If 'h', integrate data with common genes
            If 'v', integrate data profiled from the same cells
            If 'd', inetrgate data without common genes
            Default: 'h'
        loss_type
            type of loss. Choose between ['BCE', 'MSE'm 'L1']. Default: 'BCE'
        label_weight
            Prior-guided weighted vectors. Default: None
        Prior
            Prior correspondence matrices, indexed by absolute domain id: entry
            `j` is the prior between query domain `j` and the reference domain,
            so it must be long enough to be indexed at every id in
            `range(n_domain)` except `ref_id`.
        save_OT
            If True, output a global OT plan. Default: False
        use_specific
            If True, specific genes in each dataset will be considered. Default: True
        lambda_s
            Balanced parameter for specific genes. Default: 0.5
        lambda_kl: 
            Balanced parameter for KL divergence. Default: 0.5
        lambda_recon:
            Balanced parameter for reconstruction. Default: 1.0
        lambda_ot:
            Balanced parameter for OT. Default: 1.0
        reg:
            Entropy regularization parameter in OT. Default: 0.1
        reg_m:
            Unbalanced OT parameter. Larger values means more balanced OT. Default: 1.0
        lr
            Learning rate. Default: 2e-4
        max_iteration
            Max iterations for training. Training one batch_size samples is one iteration. Default: 60000
        early_stopping
            EarlyStopping class (definite in utils.py) for stoping the training if loss doesn't improve after a given patience. Default: None
        device
            'cuda' or 'cpu' for training. Default: 'cuda'
        verbose
            Verbosity, True or False. Default: False
        """

        self.to(device)

        # foreach batches the per-parameter update into grouped kernels. It is
        # the same arithmetic as the default single-tensor path (verified
        # bit-identical) but avoids one python round trip per parameter tensor,
        # which matters here because every domain adds its own decoder.
        optim = torch.optim.Adam(self.parameters(), lr=lr, weight_decay=5e-4, foreach=True)

        n_epoch = int(np.ceil(max_iteration/len(dataloader)))

        if loss_type == 'BCE':
            loss_func = nn.BCELoss()
        elif loss_type == 'MSE':
            loss_func = nn.MSELoss()
        elif loss_type == 'L1':
            loss_func = nn.L1Loss()
        else:
            raise ValueError("loss_type must be one of 'BCE', 'MSE' or 'L1', got {!r}".format(loss_type))

        # Offset of every dataset inside the concatenated index space. The
        # original code recomputed `sum(num_cell[0:j])` for every domain of
        # every mini batch.
        offsets = np.concatenate(([0], np.cumsum(num_cell))).tolist()

        query_id = list(range(self.n_domain))
        query_id.remove(self.ref_id)

        # Move the reweighting vectors to the device once instead of on every
        # mini batch (this also makes them usable with GPU-resident indices).
        if label_weight is not None:
            label_weight = [w if w is None else w.to(device) for w in label_weight]

        # Column offsets of each modality inside a vertically-concatenated cell.
        if mode == 'v':
            bounds = np.concatenate(([0], np.cumsum(num_gene[:self.n_domain]))).tolist()

        def _slice_tran(j, idx_q, idx_r):
            """Mini-batch block of the global OT plan, sliced before transfer.

            Indexing the host array first moves `len(idx_q) x len(idx_r)` values
            to the device instead of the whole `n_query x n_ref` plan.
            """
            block = tran[j][np.ix_(idx_q, idx_r)]
            return torch.from_numpy(block).to(device)

        def _slice_prior(j, idx_q, idx_r):
            """Mini-batch block of the prior, gathered without a full-width temporary.

            Callers keep the result per query domain. The previous code held a
            single `Prior_batch` variable that the domain loop overwrote, so with
            more than two domains every domain's OT used the *last* domain's
            prior block -- a shape mismatch whenever the two query domains
            contributed different cell counts, and the wrong prior when they
            happened to match.
            """
            iq = torch.as_tensor(idx_q, dtype=torch.long)
            ir = torch.as_tensor(idx_r, dtype=torch.long)
            return Prior[j][iq.unsqueeze(1), ir.unsqueeze(0)].to(device)

        with tqdm(range(n_epoch), total=n_epoch, desc='Epochs') as tq:
            for epoch in tq:

                tk0 = tqdm(enumerate(dataloader), total=len(dataloader), leave=False, desc='Iterations', disable=(not verbose))
                # Accumulate on the device in float64 and read back once per
                # epoch: calling .item() per loss per iteration forced three
                # host/device synchronisations on every step.
                epoch_loss = defaultdict(lambda: torch.zeros((), dtype=torch.float64, device=device))
                i = -1

                if mode == 'v':

                    for i, (x, idx) in tk0:
                        x = x.to(device, non_blocking=True).float()

                        x_list = [x[:, bounds[j]:bounds[j+1]] for j in range(self.n_domain)]

                        z, mu, var = self.encoder(x_list[0], 0)
                        kl_loss = kl_div(mu, var)
                        recon = self.decoder(z, 0)
                        recon_loss = loss_func(recon, x_list[0]) * 2000

                        for j in range(1, self.n_domain):

                            recon = self.decoder(z, j)
                            recon_loss += lambda_s * loss_func(recon, x_list[j]) * 2000   ## TO DO

                        loss = {'recon_loss':lambda_recon*recon_loss, 'kl_loss':lambda_kl*kl_loss}

                        optim.zero_grad(set_to_none=True)
                        sum(loss.values()).backward()
                        optim.step()

                        for k,v in loss.items():
                            epoch_loss[k] += v.detach().double()

                        if verbose:
                            tk0.set_postfix_str(','.join('{}={:.3f}'.format(k, v) for k,v in loss.items()))

                    epoch_loss = {k:(v/(i+1)).item() for k, v in epoch_loss.items()}
                    epoch_info = ','.join(['{}={:.3f}'.format(k, v) for k,v in epoch_loss.items()])
                    tq.set_postfix_str(epoch_info)


                elif mode == 'd':

                    for i, (x,y,idx) in tk0:

                        x = x.to(device, non_blocking=True).float()
                        y = y.to(device, non_blocking=True).long()
                        idx = idx.to(device, non_blocking=True)

                        if len(torch.unique(y)) < self.n_domain:
                            continue

                        mu_dict = {}
                        var_dict = {}

                        loc_ref = torch.where(y==self.ref_id)[0]
                        idx_ref = idx[loc_ref] - offsets[self.ref_id]
                        idx_ref_host = idx_ref.cpu().numpy() if (save_OT or Prior is not None) else None

                        loc_query = {}
                        idx_query = {}
                        tran_batch = {}
                        Prior_batch = {}

                        for j in query_id:

                            loc_query[j] = torch.where(y==j)[0]
                            idx_query[j] = idx[loc_query[j]] - offsets[j]

                            if save_OT:
                                idx_q_host = idx_query[j].cpu().numpy()
                                tran_batch[j] = _slice_tran(j, idx_q_host, idx_ref_host)
                            else:
                                tran_batch[j] = None

                            Prior_batch[j] = None if Prior is None else \
                                _slice_prior(j, idx_query[j].cpu().numpy(), idx_ref_host)

                        recon_loss = torch.zeros((), device=device)
                        kl_loss = torch.zeros((), device=device)
                        ot_loss = torch.zeros((), device=device)

                        loc = loc_query
                        loc[self.ref_id] = loc_ref

                        for j in range(self.n_domain):

                            # Gather the domain's rows once; the original
                            # indexed `x[loc[j]]` three times per domain.
                            x_j = x[loc[j]][:, 0:num_gene[j]]
                            z_j, mu_j, var_j = self.encoder(x_j, j)
                            mu_dict[j] = mu_j
                            var_dict[j] = var_j
                            recon_j = self.decoder(z_j, j)

                            recon_loss += loss_func(recon_j, x_j) * x.size(-1)  ## TO DO
                            kl_loss += kl_div(mu_j, var_j)

                        for j in query_id:

                            ot_loss_tmp, tran_batch[j] = unbalanced_ot(tran_batch[j], mu_dict[j], var_dict[j], \
                                mu_dict[self.ref_id].detach(), var_dict[self.ref_id].detach(), Couple=Prior_batch[j], device=device)

                            if save_OT:
                                tran[j][np.ix_(idx_query[j].cpu().numpy(), idx_ref_host)] = tran_batch[j].cpu().numpy()

                            ot_loss += ot_loss_tmp

                        loss = {'recon_loss':lambda_recon*recon_loss, 'kl_loss':lambda_kl*kl_loss, 'ot_loss':lambda_ot*ot_loss}

                        optim.zero_grad(set_to_none=True)
                        sum(loss.values()).backward()
                        optim.step()

                        for k,v in loss.items():
                            epoch_loss[k] += v.detach().double()

                        if verbose:
                            tk0.set_postfix_str(','.join('{}={:.3f}'.format(k, v) for k,v in loss.items()))

                    epoch_loss = {k:(v/(i+1)).item() for k, v in epoch_loss.items()}
                    epoch_info = ','.join(['{}={:.3f}'.format(k, v) for k,v in epoch_loss.items()])
                    tq.set_postfix_str(epoch_info)


                elif mode == 'h':

                    for i, (x, y, idx) in tk0:

                        # One host->device transfer for the whole row; the
                        # common and specific blocks are sliced on the device.
                        x = x.to(device, non_blocking=True).float()
                        y = y.to(device, non_blocking=True).long()
                        idx = idx.to(device, non_blocking=True)

                        x_c = x[:, 0:num_gene[self.n_domain]].contiguous()

                        loc_ref = torch.where(y==self.ref_id)[0]

                        idx_ref = idx[loc_ref] - offsets[self.ref_id]
                        idx_ref_host = idx_ref.cpu().numpy() if (save_OT or Prior is not None) else None

                        loc_query = {}
                        idx_query = {}
                        tran_batch = {}
                        Prior_batch = {}

                        for j in query_id:

                            loc_query[j] = torch.where(y==j)[0]
                            idx_query[j] = idx[loc_query[j]] - offsets[j]
                            tran_batch[j] = None
                            Prior_batch[j] = None

                        if len(loc_ref) > 0:
                            for j in query_id:

                                if save_OT and len(idx_query[j]) != 0:
                                    tran_batch[j] = _slice_tran(j, idx_query[j].cpu().numpy(), idx_ref_host)

                                if Prior is not None:
                                    Prior_batch[j] = _slice_prior(j, idx_query[j].cpu().numpy(), idx_ref_host)

                        ot_loss = torch.zeros((), device=device)
                        recon_loss = torch.zeros((), device=device)

                        loc = loc_query
                        loc[self.ref_id] = loc_ref
                        idx = idx_query
                        idx[self.ref_id] = idx_ref

                        z, mu, var = self.encoder(x_c, 0)
                        recon_x_c = self.decoder(z, 0, y)

                        if label_weight is None:
                            recon_loss = loss_func(recon_x_c, x_c) * 2000
                        else:
                            for j, weight in enumerate(label_weight):

                                if len(loc[j])>0:
                                    if weight is None:
                                        recon_loss += 1/self.n_domain * loss_func(recon_x_c[loc[j]], x_c[loc[j]]) * \
                                        2000
                                        # kl_loss += kl_div(mu[loc[j]], var[loc[j]])
                                    else:
                                        recon_loss += 1/self.n_domain * F.binary_cross_entropy(recon_x_c[loc[j]], x_c[loc[j]], weight=weight[idx[j]]) * \
                                        2000
                                        # kl_loss += kl_div(mu[loc[j]], var[loc[j]], weight[idx[j]])

                        kl_loss = kl_div(mu, var)
                        if use_specific:

                            x_s = x[:, num_gene[self.n_domain]:]

                            for j in range(self.n_domain):
                                if len(loc[j])>0:
                                    recon_x_s = self.decoder(z[loc[j]], j+1)
                                    recon_loss += lambda_s * loss_func(recon_x_s, x_s[loc[j]][:, 0:num_gene[j]]) * 2000

                        if len(torch.unique(y))>1 and len(loc[self.ref_id])!=0:

                            mu_dict = {}
                            var_dict = {}

                            for j in range(self.n_domain):
                                if len(loc[j])>0:
                                    mu_dict[j] = mu[loc[j]]
                                    var_dict[j] = var[loc[j]]

                            mu_ref = mu_dict[self.ref_id].detach()
                            var_ref = var_dict[self.ref_id].detach()

                            for j in query_id:
                                if len(loc[j])>0:

                                    ot_kwargs = {}
                                    if label_weight is not None:
                                        ot_kwargs['query_weight'] = label_weight[j]
                                        ot_kwargs['ref_weight'] = label_weight[self.ref_id]

                                    ot_loss_tmp, tran_batch[j] = unbalanced_ot(
                                        tran_batch[j],
                                        mu_dict[j],
                                        var_dict[j],
                                        mu_ref,
                                        var_ref,
                                        reg=reg,
                                        reg_m=reg_m,
                                        idx_q=idx_query[j],
                                        idx_r=idx_ref,
                                        Couple=Prior_batch[j],
                                        device=device,
                                        **ot_kwargs,
                                    )

                                    if save_OT:
                                        tran[j][np.ix_(idx_query[j].cpu().numpy(), idx_ref_host)] = tran_batch[j].cpu().numpy()

                                    ot_loss += ot_loss_tmp

                        loss = {'recloss':lambda_recon*recon_loss, 'klloss':lambda_kl*kl_loss, 'otloss':lambda_ot*ot_loss}

                        optim.zero_grad(set_to_none=True)
                        sum(loss.values()).backward()
                        optim.step()

                        for k,v in loss.items():
                            epoch_loss[k] += v.detach().double()

                        if verbose:
                            tk0.set_postfix_str(','.join('{}={:.2f}'.format(k, v) for k,v in loss.items()))


                    epoch_loss = {k:(v/(i+1)).item() for k, v in epoch_loss.items()}

                    epoch_info = ','.join(['{}={:.2f}'.format(k, v) for k,v in epoch_loss.items()])
                    tq.set_postfix_str(epoch_info)

                    early_stopping(sum(epoch_loss.values()), self)
                    if early_stopping.early_stop:
                        print('EarlyStopping: run {} epoch'.format(epoch+1))
                        break
