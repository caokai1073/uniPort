Installation
============

To use uniPort, first install it using pip:

.. code-block:: console

   pip3 install uniport

After a correct installation, you should be able to import the module without errors:

.. code-block:: console

   import uniport as up

We highly recommand training the model with Nvidia GPU devices, consider install Pytorch cuda version:

.. code-block:: console

   pip3 install torch torchvision torchaudio

Some of the tutorials use extra packages that uniPort itself does not depend on
(``episcanpy`` for ATAC feature selection, ``matplotlib``/``seaborn`` for the
plots). To run the notebooks, install the optional ``examples`` extra:

.. code-block:: console

   pip3 install "uniPort[examples]"