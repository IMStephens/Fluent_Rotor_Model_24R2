# Fluent_Rotor_Model_24R2
An update to the original Virtual Blade Model UDF uploaded by sfo-project on Github. The objective is to develop a rotorcraft model, with fuselage forces balanced by the main rotor, and the main rotor's torque balanced by the tail rotor.

## How To: Actuator Line Model

Disclaimer: This tutorial was written for ALM version 10.7.3, using Fluent 2024 R2.

### Geometry
As this version of the Actuator Line Model uses a floating rotor formulation, you should not designate a rotor zone in the geometry creation of your domain. It is instead advisable to create a Body of Influence in the area in which you are going to place the rotor. The dimensions of the BoI above/below, and beyond the radius of the rotor is up to the user.

### Mesh
Experiment with scaling the cell sizing in the rotor zone by the chord length of the rotor blades, as the Gaussian function that is in the ALM covers only the chord length and a smaller sponge zone. A cell sizing of 1 rotor chord is probably too coarse, so instead aim for between 3-5 cells per chord length. If using a fuselage, ensure you create a named selection for the fuselage. This will be used later in the ALM GUI. You do not need to label the BoI using named selections as the co‑ordinates of the rotor’s centre will later be defined by the user.

### Set‑up
Once the geometry and mesh are established, navigate to the working directory in which you will launch the Fluent case and ensure the correct files are present. These should be **rotor_model_vXX.c**, **rotor_model_vXX.scm**, **thread_mem.h**, **naca0012.dat** (or chosen airfoil).

Before you can begin running an ALM simulation, follow these steps within Fluent:
  - On the ribbon at the top of the Fluent application, navigate to the **User‑Defined** tab and then to **Functions**, then **Compiled…**. A window called **Compiled UDFs** should open. On the left hand side under **Source Files**, click **Add…** and find the **rotor_model_vXX.c** file and hit **OK**. Back to the **Compiled UDFs** window, on the right hand side under **Header Files** click **Add…** and find the **thread_mem.h** file. Hit **OK**. With both files appearing in the **Source Files** and **Header Files** list, tick **Use Built‑In Compiler** and then hit **Build**. The last message printed to console should be **Done**, at which point you can go back to the **Compiled UDFs** window and hit **Load**.
    
  - Again on the **User‑Defined** tab, navigate to **Function Hooks…** which should open the **User‑Defined Function Hooks** window, and next to the **Adjust** entry, press **Edit…**. You should see **my_SrcComp::libudf** appear on the **Available Adjust Functions** list. Select it and press **Add** which should move it to the **Selected Adjust Functions** list. Back to the **User‑Defined Function Hooks** window, hit **Edit…** next to the **Execute at End** entry. **my_write_rotor_data_to_file::libudf** should appear under **Available Execute at End Functions**. Select it and then click **Add** and then **OK**.
    
  - Finally, in the **User‑Defined** tab, open **Memory…** and increase the **Memory Locations** number to at least 13.
    
  - On the left‑hand side of the Fluent application, open the **Outline View**. Under **Cell Zone Conditions**, find the cell zone corresponding to the domain, which should be under **Fluid**. Double click on the cell zone to open the **Fluid** window. Ensure **Source Terms** is ticked, and then navigate to the **Source Terms** tab in the **Fluid** window. Click on **Edit…** next to **X Momentum** to open the **X Momentum sources** window and increase the **Number of X Momentum sources** to the number of rotors to simulate. There should appear a drop‑down menu for each rotor. Open the drop‑down menu next to **1.** and select **udf my_xmom_src_1::libudf**. If using a second rotor, open the drop down menu next to **2.** and select **udf my_xmom_src_2::libudf**. When finished, press **OK**. Repeat the above for **Y Momentum** and **Z Momentum**, ensuring that you select **my_ymom** for **Y Momentum**, and **my_zmom** for **Z Momentum**. Click **Apply** and **Close** on the **Fluid** window.
    
  - The Fluent simulation should now be initialised, once the user has set their own boundary conditions, turbulence models, discretization schemes etc, ensuring it is set up for a **Transient** simulation. If simulating an ambient hover, it is recommended to use **Pressure Inlet** at the inlet and **Pressure Outlet** at the outlet. These boundary conditions allow the momentum source of the rotor to drive the velocity field of the domain, but may have a longer settling time than a **Velocity Inlet** whilst the solver establishes the inlet conditions.
    
  - After initialisation, at the top of the Fluent application, navigate to **File**, then **Read**, and then **Scheme…**. Look for **rotor_model_vXX.scm** and click **OK**. Open the **Models** task on the left‑hand side of the Fluent application under **Outline View** and you should see a **Helicopter Model** entry at the bottom. Double click **Helicopter Model** to open the ALM GUI. The GUI is self-explanatory. In the **Trimming** tab, enabling **Force balance** will ignore the **Desired thrust coefficient** value in the GUI and instead calculate the vertical force on the fuselage and use this as a trim target in an attempt to balance the vertical force generated by the rotor with the fuselage drag, and weight. Enabling **Moment balance** will use the torque produced by the main rotor (rotor 1) as a thrust target for the tail rotor (rotor 2). The code will find the distance between the 2 rotors and divide the main rotor torque by that distance to find the required thrust value. The **Define CG Position** is currently defunct. When you have defined a rotor, hit **Change/Create** before changing the value of **Active Rotor Zone** to the next rotor. Once all rotors are defined, click **OK**.
    
  - At this stage, you are ready to begin iterating. 

## Cautions
There are a few things to consider when setting up and running an ALM solution. 
  - Chord and twist are linearly interpolated between sections.
  - There is a limit of 20 sections per rotor zone.
  - ALM only supports main rotor or main & tail rotor configurations. It is not currently configured for quadcopters for example, but this could be built.
  - Rotor 1 must be the main rotor and rotor 2 the tail rotor.
  - The coefficient values are calculated in the EU form, i.e: T = 0.5 * rho * (V²) * S * C_T
