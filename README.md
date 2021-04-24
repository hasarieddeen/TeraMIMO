# ---- TeraMIMO Channel Simulator v1.0. ----


 -- A Channel Simulator for Wideband Ultra-Massive MIMO Terahertz Communications


 -- (c) 2021 Simon Tarboush, Hadi Sarieddeen, Hui Chen, Mohamed Habib Loukil, Hakim Jemaa, 
             Mohamed-Slim Alouini, Tareq Y. Al-Naffouri

 -- e-mail: simon.w.tarboush@gmail.com; hadi.sarieddeen@kaust.edu.sa; hui.chen@kaust.edu.sa;
            mohamedhabib.loukil@kaust.edu.sa; hakim.jemaa@kaust.edu.sa;
            slim.alouini@kaust.edu.sa; tareq.alnaffouri@kaust.edu.sa




----  Important information ----


If you are using this simulator (or parts of it) for a publication, then you must cite the following paper:

S. Tarboush, H. Sarieddeen, H. Chen, M.-H. Loukil, H. Jemaa, M.-S. Alouini, and T. Y. Al-Naffouri, "TeraMIMO: A channel simulator for wideband ultra-massive MIMO terahertz communications," arXivpreprint arXiv:2104.11054, 2021


For more on signal processing for THz communications applications that could utilize TeraMIMO, please check the paper:

H. Sarieddeen, M.-S. Alouini, and T. Y. Al-Naffouri, "An overview of signal processing techniques for terahertz communications," arXivpreprint arXiv:2005.13176, 2020





----  How to start a TeraMIMO simulation ----


Simply run Generate_TIV_THz_Channel.m (or Generate_TV_THz_Channel.m for the time-variant version)

Configure the required simulation parameters in generate_channel_param_TIV.m (such as channel type,
operating frequency and bandwidth, AoSA parameters, geometry, and THz-specific features of the channel, the spherical wave model, beam split, and misalignment)

The code computes the molecular absorption coefficient using compute_Abs_Coef.m based on the selected method in generate_channel_param_TIV.m

The main channel computations for both frequency and delay domains are executed in channel_TIV.m

For more details on the inputs and outputs, check the "help" text of each function

You can visualize the channel using Plot_TIV_THz_Channel

We highly recommend you execute the code step-by-step (using MATLAB's debug mode) to get a detailed understanding of the simulator





----  How to run GUI ----


1- Go to folder "GUI" and open TeraMIMO_GUI.m; run it in MATLAB

2- In the GUI, select the parameters as required for your simulation and click on “Run"!

3- The results are saved to the workspace!




Note: The code assumes the presence of a few gases in the medium in its default mode. You can experiment with more combinations of gases by adding more .csv files for different gases from the folder "Molecular_Absorption," subfolder "Data_Full," to the subfolder "Data,". Due to the size limitation on GitHub, "Data_Full" does not include all HITRAN gases. Feel free to reach out to send you the rest.



----  Version history ----

* Version 1-0: initial version for GitHub release
