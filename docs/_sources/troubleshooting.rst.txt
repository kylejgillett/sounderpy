❓ Troubleshooting
===================


Like any software, SounderPy is prone to errors now and then -- some of them are known with fixes! If you are running into issues using SounderPy, check out this guide to troubleshooting below! 

If you have an error that you don't see on this page, you can `open a GitHub Issue <https://github.com/kylejgillett/sounderpy/issues>`_ or feel free to shoot me a DM on `Twitter <https://twitter.com/wxkylegillett>`_!

Funky Plot Scaling 
-------------------

.. image:: _static/images/funky_layout_example.jpg
   :width: 300 px
   :align: left

When using SounderPy's plot function, sometimes the matplotlib formatting operations the plot functions use (``build_sounding()``, ``build_hodograph()``, ``build_composite()``) don't agree with some screen sizes. This throws off the scaling of the plots when plotting them 'inline'. Thankfully there is a fix! 

Each of SounderPy's plotting functions have a ``save`` key-word-argument (kwarg). By default these are set to *False*, but when setting the kwarg to *True* you can save the plot to a .png file. 

When the plot is saved to a file it will preserve the correct scaling and dimensions. 

The Fix: saving the plot to a file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	Add the ``save=True`` kwarg to your plot function (``build_sounding()``, ``build_hodograph()``, ``build_composite()``). 

	You can optionally add the ``filename`` kwarg that allows you to set a specific filename. By default, this filename is *'sounderpy_sounding.png'*.