import numpy as np

def scan_timesteps(file):
    timesteps = []
    
    with open(file, 'r') as f:
        text = f.read()
    lines = text.splitlines()
    
    
    for l_num, line in enumerate(lines):
        line = line.split()
        if (len(line) ==2) and (line[1] == "TIMESTEP"):
            timesteps.append(lines[l_num+1])
    return timesteps

def read_dump(file, ret_dict = False, tsi = None, atomic = False, boxdims = False, skips_ixs = False, atomic_forces = False):
   
        
    with open(file, 'r') as f:
        text = f.read()
    lines = text.splitlines()
    
    ts_scan = scan_timesteps(file)
    if np.any(tsi):
        ts_incl = np.array(ts_scan)[np.array(tsi)]
    else:
        ts_incl = np.array(ts_scan) 
    save_timestep = False
    
    atoms = int(lines[3].split()[0])
    ts = 10000 #if you every have a trajectory with more than this many timesteps, this will need to be changed
    xpos = np.zeros([ts, atoms])
    ypos = np.zeros([ts, atoms])
    zpos = np.zeros([ts, atoms])
    
    timesteps = []
    len_box = []
    atype= np.zeros([atoms], dtype = 'int')
    molecnum= np.zeros([atoms], dtype = 'int')
    
    count = 0
    
    if atomic:
        ndata = 5
        shift = 0
    elif atomic_forces:
        ndata = 8
        shift = 3
    else:
        ndata = 6
        shift = 0
        
    tscount = -1
    for l_num, line in enumerate(lines):
        line = line.split()
        if (len(line) ==2) and (line[1] == "TIMESTEP"):
            
            ts = lines[l_num+1]
            tscount+=1
            if ts in ts_incl:
                save_timestep = True
            else:
                save_timestep = False
            if save_timestep:
                timesteps.append(ts)
            
        elif (len(line) ==6) and (line[1] =="BOX"):
            a = lines[l_num+1].split()[1]
            if save_timestep:
            	if boxdims:
            		len_box.append([float(x) for x in lines[l_num+1].split()])
            	else:
                	len_box.append(float(a))
                
        elif (len(line) ==ndata) and (line[0] != "ITEM:"):  
            if save_timestep:
                if skips_ixs:
                    row = int(count/atoms)
                    col = (count - atoms*int(count/atoms))
                else:
                    a = int(line[0])
                    row = tscount
                    col = a - 1

                    
                x = float(line[ndata-3-shift])
                y = float(line[ndata-2-shift])
                z = float(line[ndata-1-shift])            
                if count<atoms:
                    atype[count]= int(line[ndata-4-shift])
                    if not atomic and not atomic_forces:
                        molecnum[count] = int(line[1])
                
                # xpos[int(count/atoms), (count - atoms*int(count/atoms))]=x
                # ypos[int(count/atoms), (count - atoms*int(count/atoms))]=y
                # zpos[int(count/atoms), (count - atoms*int(count/atoms))]=z
                xpos[row, col]=x
                ypos[row, col]=y
                zpos[row, col]=z
                
                count += 1
            
        else:
            continue
    
    xcnt = 0
    ycnt = 0
    zcnt = 0
    if np.all(xpos[0,:]==0):
        xcnt =1
    if np.all(ypos[0,:]==0):
        ycnt = 1
    if np.all(zpos[0,:]==0):
        zcnt = 1
        
    xpos = xpos[~np.all(xpos == 0, axis=1)]
    ypos = ypos[~np.all(ypos == 0, axis=1)]
    zpos = zpos[~np.all(zpos == 0, axis=1)]
    timesteps = [int(ts) for ts in timesteps]
    if xcnt ==1:
        xpos = np.vstack((np.zeros(atoms), xpos))
    if ycnt ==1:
        ypos = np.vstack((np.zeros(atoms), ypos))
    if zcnt ==1:
        zpos = np.vstack((np.zeros(atoms), zpos))
        
    if ret_dict:
        ret = {'box_lens': len_box,
               'timesteps': timesteps,
               'xpos': xpos,
               'ypos':ypos,
               'zpos':zpos,
               'types':atype,
               'molecs':molecnum,
               'n_atoms':atoms
            }
        return ret
        
        
        
    return (len_box, timesteps, xpos, ypos, zpos, atype, molecnum, atoms)
