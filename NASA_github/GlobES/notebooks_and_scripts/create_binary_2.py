import numpy as np
import os
import time
        
folder= '/Users/apayne3/Desktop/GlobES/Checlair_plots/'


# os.listdir() method in python is used to get the list of all files and directories in the specified directory. If we don't specify any directory, then list of files and directories in the current working directory will be returned
for file in sorted(os.listdir()):
    if 'atm.txt' in file:
        open_file= open(folder+file, 'r')
        atm_file= open_file.readlines()
        # close the file!
        open_file.close()

        # GOAL: create arrays (for the 3D model) filled with atmospheric data extracted from the atmospheric layering files

        # Gather the names of the molecules contained in each atmosphere file
        molecules=[]
        for line in atm_file:
            if '<ATMOSPHERE-LAYERS-MOLECULES>' in line:
                molecule= line.split(',')
                for value in molecule:
                    molecules.append(value)
                    # Reformatting to get rid of the excess data contained in the slicing
        molecules[0]= molecules[0][-3:]
        if 'E0' in file:
            molecules[-1]=molecules[-1][-5:-1]
        else:
            molecules[-1]= molecules[-1][-6:-1]

        # Here we gather all of the abundance information from ATMOSPHERE-LAYER- and store it in one long list  
        test=[]
        for line in atm_file:
            if "<ATMOSPHERE-LAYER-" in line:
                # Add 2 to the # of entries bc the first two values in these layers are pressure and temp
                # which are necessary for GlobES
                values= np.arange(0,len(molecules)+2,1)
                for value in values:
                    result= line.split(',')[value]
                    if '>' in result:
                        result= result.split('>')[1]
                    test.append(float(result)) 
        

        # store the # of layers in each file (as integer)
        for line in atm_file:
            if "<ATMOSPHERE-LAYERS>" in line:
                num_layers= line.split('>')[1]
                num_layers_int= int(num_layers)

        # extract the atm abund based on its index value and REORDER them...
        # This puts them in a sequence where the lists of each molecule in MOLECULES will comes one after the next
        layers= np.arange(0,num_layers_int,1)
        data_lists=[]
        for value in values:
                sorted_data=[]
                for layer in layers:
                    # add 2 because the values in test now include temp and pressure
                    result=test[value+((len(molecules)+2)*layer)]
                    sorted_data.append(result)
                data_lists.append(sorted_data)

        # This separates the long list that is data_lists and creates a bunch of smaller lists of the abundance values in order
        # this way you can index an abundance_file and all of the values will correspond to a atm abundance list of a given molecule (and the first two are Press and Temp)
        abundance_files=[]
        for data in data_lists:
            abund_chunk=[]
            for i in np.arange(0, len(data), (num_layers_int+2)):
                abun_chunk= data[i:i + (2+num_layers_int)]
                abund_chunk.append(abun_chunk)
            abundance_files.append(abund_chunk)
        # print(np.shape(abundance_files))
        # print(abundance_files[2][0])
        
        # including the start and end time only for reference to improve speed of the for loop that creates an array
        start = time.time()
        # 28 times
        # creates an array for: 1 of each input files, all (28) molecules, one for each layer, 144 lat points, 91 long points
         # it should create 6x28 resultant arrays with 144x91 data points
         # where every point in array is equal to the value of the molecule at layer that layer
         # beginning with a zeros array and replacing the values as you go is a much faster process than looping over every value in array 
        store_data_vk = np.zeros(shape=(len(molecules),num_layers_int,144,91))
        for i in range(len(molecules)):
            rep_arr = np.repeat(abundance_files[i][0],repeats=144*91)
            rep_arr_rs = np.reshape(rep_arr,newshape=(num_layers_int,144,91))
            store_data_vk[i,:,:,:] = rep_arr_rs = np.reshape(rep_arr,newshape=(num_layers_int,144,91)) #np.tile(abundance_files[i][0],reps=(num_layers_int,144,91))

        #print(store_data_vk[0])
        
        end = time.time()
        print("The loop took" , end - start)
        
        # Export to binary format here...
        file_list=[]
        if 'E0' in file:
            name= file[0:6]
            file_list.append(name)
        else:
            name= file[0:7]
            file_list.append(name)
        
        for file in file_list:
            outfile = 'binary_file_out_'+file+'.dat'
            # the code will raise an error if the file already exists- it wont add data to it
            # this line will delete the file and create a new one if it already exists
            try: os.remove(outfile)
            except: pass
            with open(outfile, 'xb') as fo:
                #for entry in store_data:
                fo.write(bytes('<BINARY>',encoding = 'utf-8'))
                fo.write(np.array(store_data_vk,order='C'))
                fo.write(bytes('</BINARY>',encoding = 'utf-8'))

        break