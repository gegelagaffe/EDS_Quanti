#Function to provide interresting information for SEM-EDS analisys....V0.1

def Carac_Mineral(Mineral, Elt, Mol, rho, E0, Opt, Graph):
    
    import matplotlib.pyplot as plt
    import hyperspy.api as hs
    import numpy as np
    
    if Mineral == 'Help':
        Output0 = Mineral
        
        Mineral = 'Portlandite'
        Elt_CH= ['H', 'Ca', 'O']
        Mol_CH= [2,1,1]
        rho_CH= 2.24
        E0= 15
        Out_Elt = 'Tot'
        Graph= 'Y'
        
        print('Function which provides interresting informations for SEM-EDS analisys ')
        print('')
        print('   A22 = Carac_Mineral(Nom_CH, Elt_CH, Mol_CH, rho_CH, E0, Out_Elt, Graph)')
        print('')
        print('Name of the mineral:.......................................  ','Nom_CH =', Mineral)
        print('Ionic Species:.............................................  ','Elt_CH =', Elt_CH)
        print('Numbre of mol of Each Ionic Species (mol):.................. ','Mol_CH =', Mol_CH)
        print('Density of the mineral (g/cm3):............................. ','rho_CH =', rho_CH)
        print('Energy of SEM analysis (KeV):..............................  ','E0 =', E0)
        print('Output otpion (Tot: Global information Ca : only calcium) :  ','Out_Elt =', Out_Elt)
        print('Print option Y:print otherwise no print:...................  ','Graph =', Graph)        
        
    else:
        #INITIALISATION
        ###############
        #Number of ionic species
        dim_IS = len(Elt)
        #Name of ionic species
        Name_IS = list(range(dim_IS))
        #atomic weight
        A_IS = list(range(dim_IS))
        #Atomic number
        Z_IS = list(range(dim_IS))
        #Xray_line_Ka
        XR_Ka = list(range(dim_IS))
        #Coefficient d'absorption massique µ/rho cm^2/g  
        muRO_Ka = list(range(dim_IS))
        muRO_Ka[0] = 0
        #Properties ionic species
        PropIS = list(range(dim_IS))
        #General properties of the mineral
        #density (g/cm3)
        rho_IS = rho
        # Calcul of depth penetration of X ray (disturbance peer)
        SVAR0 = 0
        #Tension d'acceleration (KeV)
        EtA = E0

        #CALCULATION
        ###############
        #massic ratio
        C_Elt = hs.material.atomic_to_weight(Mol, Elt)
        C_Elt = C_Elt.tolist()
        # Atomic ratio
        N_Elt = hs.material.weight_to_atomic(C_Elt, Elt)
        N_Elt = N_Elt.tolist()
        #
        #For mass absorption calculation
        position_H = Elt.index('H')
        Elt2 = Elt[:]
        C_Elt2 = C_Elt[:]
        del Elt2[position_H]
        del C_Elt2[position_H]
        #
        #
        #For loop to calculate name, Z, Atomic weight Ka_Xrayline and mass absorbtion coefficient
        for i in range(dim_IS):
            PropIS[i] = hs.material.elements[Elt[i]].General_properties 
            Name_IS[i] = PropIS[i].name
            Z_IS[i] = PropIS[i].Z
            A_IS[i] = PropIS[i].atomic_weight
            XR_Ka[i] = hs.material.elements[Elt[i]].Atomic_properties.Xray_lines.Ka['energy (keV)']
            VAR0 = C_Elt[i]*0.01*A_IS[i]/Z_IS[i]
            SVAR0+=VAR0
            if i >  position_H:
                muRO_Ka[i] = hs.material.mass_absorption_mixture(elements= Elt2, weight_percent= C_Elt2, energies=[ XR_Ka[i]])[0]
        #
        #
        Depth= 0.033*((EtA**(1.7)) / rho_IS)*SVAR0
        Depht_2 = [Depth]*dim_IS
        #
        #output : if you propose Total, the output would be the whole matrix whereas if you select just one element you'll get the elt informations
        if Opt=='Tot':
            Output0 = [Mineral, Name_IS, Z_IS ,N_Elt , C_Elt, XR_Ka, muRO_Ka, Depth ]
        else:
            if Opt in Elt:
                pos_Opt = Elt.index(Opt)
                Output0 = [Mineral, Name_IS[pos_Opt],Z_IS[pos_Opt],N_Elt[pos_Opt],C_Elt[pos_Opt], XR_Ka[pos_Opt], muRO_Ka[pos_Opt], Depth ]
            else:
                print('The element choosen is not present in the mineral')
                Output0 = [Mineral, 'Error','Error','Error','Error', 'Error', 'Error', 'Error' ]
        if Graph == 'Y':
            print('Name of Mineral.....................',Output0[0])
            print('Name of Elt.........................',Output0[1])
            print('Atomic number.......................',Output0[2])
            print('Atomic ratio in the mineral.........',Output0[3], '%')
            print('Massic ratio in the mineral.........',np.around(Output0[4],3), '%')
            print('Energy of Xray line Ka..............',np.around(Output0[5],3), 'KeV')
            print('Mass absorption coefficient µ/rho...',np.around(Output0[6],3), ' cm^2/g')
            print('Disturbance peer....................',round(Output0[7],2), 'µm')
    
    return Output0