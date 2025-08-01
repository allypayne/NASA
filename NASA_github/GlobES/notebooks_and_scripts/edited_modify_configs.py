import os
import numpy as np

# Path to folder containing the 5 atmosphere files from GV
# Files E-0 to E-5
folder='/Users/allypayne/Desktop/NASA_GSFC/IntroPlots/FilesfromGV/'
# direction to the file
baseconfig= open("/Users/allypayne/Desktop/NASA_GSFC/ConfigCode/psg_cfg_140.txt", "r");
# reads/ stores each line of file in basecfg (as a list)
basecfg= baseconfig.readlines();
#close the file everytime
baseconfig.close();

## (Here we find "c1"- line where the base_config needs to be changed by the atm files)
# BASE CONFIG: "<ATMOSPHERE-LAYERS-MOLECULES>H2O,CO2,O3..."
c_base=0
for line in basecfg:
    if "ATMOSPHERE-LAYERS-MOLECULE" in line:
        #print(c_base, line, "c_base")   
        break
    c_base= c_base+1
# to fix the python indexing    
#c_base=c_base+1

# BASE CONFIG: "<ATMOSPHERE-ABUN>1,1,1,1,1,1,1,1"
loc_abun_cfg=0
# Create a counter to find location of the abundance values in the config template
for line in basecfg:
    if "ATMOSPHERE-ABUN" in line:
        #print+ break only execute when the if statement is satisfied
        #print(loc_abun_cfg, line, "loc_abun_cfg")
        break
    loc_abun_cfg= loc_abun_cfg+1
loc_abun_cfg=loc_abun_cfg+1

# BASE CONFIG: "<ATMOSPHERE-LAYERS>60"
# Create a counter to find location of the # layers line in the base config file
loc_num_layers_base=0
for line in basecfg:
    if "<ATMOSPHERE-LAYERS>" in line:
        #print(loc_num_layers_base, line, "loc_num_layers_base")
        break
    loc_num_layers_base= loc_num_layers_base+1
#fix python indexing
#loc_num_layers_base= loc_num_layers_base+1
print(loc_num_layers_base)

# BASE CONFIG: location where "<ATMOSPHERE-LAYER-1>" is mentioned
loc_atm_layers_base=0
for line in basecfg:
    if "<ATMOSPHERE-LAYER-1>" in line:
        #print(loc_atm_layers, line, "loc_atm_layers")
        break
    loc_atm_layers_base= loc_atm_layers_base+1
#loc_atm_layers_base=loc_atm_layers_base+1
print(loc_atm_layers_base)

# INFO from E- Files from Geronimo
#list of strings for each input file name
list_atm= ['_all', '_no_o2', '_no_o3']
for file in sorted(os.listdir(folder)):
    #ignore the DS files
    if '.DS' in file: continue
    # open all of the env files from geronimo
    e_file= open(folder+file, "r");
    elines= e_file.readlines();
    # close the e_file here and continue to use elines for the rest of the code...
    e_file.close()

    # E-File: "<ATMOSPHERE-LAYERS>60"
    # Create a "counter" for the eline (will reset with every atmosphere file thats opened)
    loc_num_lines_efile=0
    for line in elines:
        # stores the # of layers in the file as a string
        if r'<ATMOSPHERE-LAYERS>' in line:
            # find the # of layers from the atm file and convert it from a string to integer
            num_layers= int(float(line.split('>')[-1]))
            print(loc_num_lines_efile, line, "loc_num_lines_efile", num_layers)
    loc_num_lines_efile=loc_num_lines_efile+1
    print(num_layers)    

    # E-File: "<ATMOSPHERE-LAYERS-MOLECULES>H2O,CO2,O3..."
    c_e_molecules=0
    for line in elines:
        # find the molecules that we include in the atmosphere from e file
        if "ATMOSPHERE-LAYERS-MOLECULE" in line:
            print(c_e_molecules, line, "c_e_molecules")
            new_molecules= line
            #include break so it stops once it finds it
            break
        c_e_molecules= c_e_molecules+1
    c_e_molecules= c_e_molecules+1    
    
     # E-File: "<ATMOSPHERE-ABUN>1,1,1,1,1,1,1,1"
    loc_e_abun= 0    
    for line in elines:
        if "ATMOSPHERE-ABUN" in line:
            print(loc_e_abun+1, line, "loc_e_abun")
            break
        loc_e_abun= loc_e_abun+1
    loc_e_abun=loc_e_abun+1

    loc_atm_layers_efile=0
    for line in file:
        if "<ATMOSPHERE-LAYER-1>" in line:
            #print(loc_atm_layers_efile, line, "loc_atm_layers_efile")
            break
        loc_atm_layers_efile= loc_atm_layers_efile+1
        #unsure why i need to add 2 here...
    loc_atm_layers_efile= loc_atm_layers_efile+1

    #make all modifications to the template file here
    #this will be executed 6 times (once for each file)
    for atm in list_atm:
        print(atm)
        # this will be executed 3 times (total of 18 files)
        file_name= '/Users/allypayne/Desktop/NASA_GSFC/IntroPlots/NewPSGFiles_140/'+file[0:7]+atm
        # open file in write mode to store the new config files
        f = open(file_name, 'w')
        # give each atm a local config file to modify (use the read lines version)
        new_file= basecfg.copy()
        print(file_name)

        ## adding this to see what the file looks like
        # CHECKED here and the new file contains all of the info (including lmax, nmax)
        # for line in new_file:
        #     test= f.write(line)
        # f.close()

        # This replaces the line of list of molecules with line with the unique atm file line
        new_file[c_base]= new_molecules
        # all checked to here

        # Change abundance values for all 3 types:
        new_file[loc_num_layers_base]= "<ATMOSPHERE-LAYERS>"+str(num_layers)+"\n"

        for i in range(num_layers):
            try: new_file[loc_atm_layers_base+i]= elines[loc_atm_layers_efile+i]
                # this will overwrite the data that exists in the file already until it has overwritten everything there
            except: new_file.append(elines[loc_atm_layers_efile+i])
                # if there is more data in the new field it will append the rest at the end

#list_atm= ['_all', '_no_o2', '_no_o3']
# THESE LINES DONT CHANGE THE INPUTS
        # no o2
        loc_abun_cfg=0
# Create a counter to find location of the abundance values in the config template
        for line in basecfg:
            if "ATMOSPHERE-ABUN" in line:
                #print+ break only execute when the if statement is satisfied
                #print(loc_abun_cfg, line, "loc_abun_cfg")
                break
            loc_abun_cfg= loc_abun_cfg+1
        #loc_abun_cfg=loc_abun_cfg+1 

        for line in new_file: 
            f.write(line)
        f.close()

        # ALL
        if atm  == list_atm[0]:
            #for line in new_file: 
            # The config file only uses the last entry for abundance value so we can just append
            # using geronimos appending code for abundance values:
            nf= open(file_name, 'a') ; nf.write('<ATMOSPHERE-ABUN>1,1,1,1,1,1,1,1'); nf.close()
        print("code passed first if statement")
        
        # NO O2
        if atm == list_atm[1]:
            #for line in new_file:
            print('changing abundance line to atmosphere', atm)
            #new_file[loc_abun_cfg]== '<ATMOSPHERE-ABUN>1,1,1,1,1,1,0,1'
            nf= open(file_name, 'a') ; nf.write('<ATMOSPHERE-ABUN>1,1,1,1,1,1,0,1'); nf.close()
        print("code passed second if statement")

        #NO O3
        if atm == list_atm[2]:
            #for line in new_file:
            nf= open(file_name, 'a') ; nf.write('<ATMOSPHERE-ABUN>1,1,0,1,1,1,1,1'); nf.close()
            #print(new_file[loc_abun_cfg]) 
            #    break

        print(atm)
        
        # for line in new_file: 
        #     f.write(line)
        # f.close()

       # print(loc_abun_cfg)

        # for line in new_file: 
        #     f.write(line)
        # f.close()


## The code below will give the string file needed to call the API and run the files created directly by running this file
folder_2= '/Users/allypayne/Desktop/NASA_GSFC/IntroPlots/NewPSGFiles_140/'
#run_this=''
for file in sorted(os.listdir(folder_2)):
    if '.ipynb' in file: continue
    if 'DS' not in file:
        file_str= str(file)
        value= 'curl --data-urlencode file@'+folder_2+file_str+' https://psg.gsfc.nasa.gov/api.php > ./IntroPlots/PSG_140_output/' +file_str+'_out'
        # use this is .append but for strings instead of list objects
        os.system(value)
        # using break only tests the first one
        #break
        #run_this += value

