## A program which reads a text file with lattice information
## and then create a stim quantum circuit.
## written by Renyu Wang (UCR)
import numpy as np
import re
from collections import Counter
# import stim (not actually needed in this function)

def do_code(filename:str, p1:float, p2:float): #fullname with extension e.g.: ‘readme.txt’
    ## default error rates
    p = [p1, p2, 0, 0, 0] #10^-3 - 10^-2
    #p = [0, 0, 0, 0, 0]

    ## read file
    na = 0
    nl = 0
    nd = 0
    array = []
    parray = []
    matrix = []
    pmatrix = []
    ls = ""
    savelnums = []
    saved = []
    dnum = []
    dmatrix = []
    
    with open(filename, 'r') as f:
        alltxt = f.readlines()    # read all text
        for line in alltxt:       # for each line
            split = line.split()  # split line by space
            if len(split) != 0: 
                if split[0] != "#":   
                    if split[0].isnumeric() == True:   # if the first character is a number, this line writes the size dr*dc 
                        ndat = int(split[0]) # how many data qubits
                        rep = int(split[1])  # how many measurement repeations
                        maxop = int(split[2]) # maximum operations on one ancilla qubit
                        #print("# n =",ndat,", rep =",rep,", maxop =",maxop,"\n")
                        
                    elif split[0] == "!": #if the first character is "!", this is for decetor evernts
                        for i in range(1,len(split)):
                            if split[i] == "#":    # if there is a #, end the line
                                break
                            else:
                                dnum.append(int(split[i]))
                        nd += 1          # add a detector event
                        dmatrix.append(dnum)
                        dnum = []
                    elif split[0] == ":":  # if the first character is ":", this is a logical operator 
                        for i in range(len(split)):
                            if split[i] == "#":    # if there is a #, end the line
                                break
                            else:
                                ls += split[i]           # add the whole line
                        xy = re.findall(r"X|Y",ls)       # find if there is "X" or "Y" exist
                        if len(xy) == 0:                 # if all Z
                            nl += 1                      # this logical operator is a detector event
                            lnums = re.findall(r"\d+", ls)     # find all numbers in the line
                            savelnums.append(lnums)            # save numbers for detector event
                        ls = ""
                        
                              
                    else:                 # the rests are stablizers
                        na += 1           # each stablizer needs an ancilla qubit
                        for i in range(len(split)):
                            #print("i=",i,"split[i]=",split[i],"\n")
                            if split[i] == "#":    # if there is a #, end the line
                                break
                            elif i < maxop:        # the item describes the stablizer
                                array.append(split[i])
                            else:                  # the item describes the error rate
                                parray.append(ord(split[i])-65)
                     
                        if len(parray) != 0:      # record if there is a special error rate 
                            pmatrix.append(parray)
                        else:
                            pmatrix.append([0]*maxop)
                    
                        matrix.append(array)  
                        array = []
                        parray = []

    Matrix = np.array(matrix)                    # transpose 
    Pmatrix = np.array(pmatrix)
    Tmatrix = Matrix.T
    Tpmatrix = Pmatrix.T
    

    ## check if there are two gates acting simultaneously on a qubit
    s = ""
    for i in range(maxop):
        for j in range(na):
            s += Tmatrix[i][j]           # add the whole line
        nums = re.findall(r"\d+", s)     # find all numbers in the line
        s = ""
        if len(set(nums))!=len(nums):    # if there are duplicates 
            print(f"two gates acting simultaneously on a qubit at step {i+1}")
            exit()

    ntot = ndat + na
    
    
    ## get last round detector events
    ds = ""
    dlast = []
    dlastmatrix = []
    clist = []
    cmatrix = []
    for i in range(nd):
        for j in range(len(dmatrix[i])):
            for k in range(maxop):
                clist.append(Matrix[dmatrix[i][j]-1][k])
                ds += Matrix[dmatrix[i][j]-1][k]
            dnums = re.findall(r"\d+", ds)
            for ele in dnums:
                dlast.append(int(ele))
            dlast.append(-na-1+dmatrix[i][j])
            ds = ""
        ### prepare for commutivite check
        #print(i,"clist=",clist)
        Clist = []
        for l in range(len(clist)):
            if clist[l] not in Clist:
                Clist.append(clist[l])
            else:
                Clist.remove(clist[l])
        #print(i,"Clist=",Clist)
        cmatrix.append(Clist)
        clist = []
        ### prepare for last round detector events
        #print(i,"dlast=",dlast)
        Dlast = []
        for m in range(len(dlast)):
            if dlast[m] not in Dlast:
                Dlast.append(dlast[m])
            else:
                Dlast.remove(dlast[m])
        #print(i,"Dlast=",Dlast)
        dlastmatrix.append(Dlast)
        dlast = []
    #print(dlastmatrix)
    
    ## check if the detect conserved generators commute with all normal generators
    string = ""
    for k in range(len(cmatrix)):
        for l in range(len(cmatrix[k])):
            string += cmatrix[k][l]+"-"
        for i in range(na):
            check = 0
            for j in range(maxop):
                if "Z" not in matrix[i][j] and "_" not in matrix[i][j]:
                    cnum = re.findall(r"\d+",matrix[i][j])
                    cnumz = "Z%s-"%(cnum[0])
                    if cnumz in string:
                        check += 1
            if check%2 != 0:
                print(f"g{dmatrix[k]} dose not commute with g{i+1}" )
                exit()
        string = ""
        
        
    ## separete X Y Z gates and get the numbers of operation qubits         
    x = "" #XCX from data to ancilla
    z = "" #CX from data to ancilla
    y = "" #YCX from data to ancilla
    nline = []
    nmatrix = []
    allx = []
    allz = []
    ally = []
 
    for i in range(maxop):
        for j in range(na):
            if "X" in Tmatrix[i][j]:         
                num = re.findall(r"\d+", Tmatrix[i][j])
                for ele in num:
                    x += ele+" "+str(ndat+j)+" "        # add numbers of qubits with X operation 
                    nline.append(str(ndat+j)+" "+ele+" ")
            elif "Z" in Tmatrix[i][j]:
                num = re.findall(r"\d+", Tmatrix[i][j])
                for ele in num:
                    z += ele+" "+str(ndat+j)+" "        # add numbers of qubits with Z operation
                    nline.append(ele+" "+str(ndat+j)+" ")
            elif "Y" in Tmatrix[i][j]:
                num = re.findall(r"\d+", Tmatrix[i][j])
                for ele in num:
                    y += ele+" "+str(ndat+j)+" "        # add numbers of qubits with Y operation
                    nline.append(ele+" "+str(ndat+j)+" ")
            else:
                nline.append("")
                
        nmatrix.append(nline)                 # save numbers of qubits matrix for each operation
        nline = []
        allx.append(x)
        x = ""
        allz.append(z)
        z = ""
        ally.append(y)
        y = ""
    #print("XCX:",allx,"\n CX:",allz,"\n YCX:",ally,"\n")
    #print(nmatrix,"\n")
        
    ## get error rates
    def get_error(ind: str):
        e = ""
        es = []
        if p[0]>0:
            e0 = f"\n{ind}DEPOLARIZE2({p[0]}) "
        else:
            e0 = ""
    
        for i in range(maxop):
            for j in range(na):
                if Tpmatrix[i][j] == 0:      # if p[0] (most cases)
                    if p[0]>0:
                        e0 += nmatrix[i][j]                     
                else:                        # else, write one by one
                    if p[Tpmatrix[i][j]]>0:
                        e += f"\n{ind}DEPOLARIZE2({p[Tpmatrix[i][j]]}) "+nmatrix[i][j]
            es.append(e0+e)
            if p[0]>0:
                e0 = f"\n{ind}DEPOLARIZE2({p[0]}) "
            else:
                e0 = ""
            e = ""
        return es
    
    
    ## begin the actual construction of the code
    sout="# begin circuit\n"
    indent="    "
    
    # beg of def do_round() "{"
    def do_round(sout: str):
        """Insert one full measurement round.
        """
        sout +=f"\n{indent}TICK"
    
        #for i in range(maxop):            # for each step
            #if allx[i] != "":             # write X Z Y gates respectively
                #sout += f"\n{indent}XCX "
                #sout += f"{allx[i]}"
            #if allz[i] != "":
                #sout += f"\n{indent}CX "
                #sout += f"{allz[i]}"
            #if ally[i] != "":
                #sout += f"\n{indent}YCX "
                #sout += f"{ally[i]}"
            #error = get_error(indent)
            #sout += error[i]                 #insert error
        
        #for subsystem code: separate X and Z measurements
        sout += f"\n{indent}XCX "
        for i in range(len(allx)):
            sout += f"{allx[i]}"
        
        sout += f"\n{indent}CX "
        for i in range(len(allz)):
            sout += f"{allz[i]}"
        
        error = get_error(indent)
        for i in range(maxop):
            sout += error[i]                 #insert error
           
        sout +=f"\n{indent}TICK"
        
        return sout
        # "}" end of def do_round()
 
    
    # beginning of circuit
    # reset all qubits to [0] state
    sout += f"{indent}R"               
    for i0 in range(ntot):
        sout += f" {i0}"

    if p[0]>0:
        sout += f"\n{indent}DEPOLARIZE1({p[0]}) "
        for i0 in range(ntot):
            sout += f" {i0}"

    # insert a measurement round
    sout=do_round(sout)
    
    if p[0]>0:                        # depolarize just before measurement
        sout += f"\n{indent}DEPOLARIZE1({p[0]}) "
        for i0 in range(ndat,ntot):
            sout += f" {i0}"

    sout += f"\n{indent}MR"
    for i0 in range(ndat,ntot):
        sout += f" {i0}"
    
    # insert 1st round detector events
    hnd = int(nd/2)
    for i in range(hnd,nd):
        sout += f"\n{indent}DETECTOR"
        for j in range(len(dmatrix[i])):
            sout +=f" rec[{-na-1+dmatrix[i][j]}]"
             
    # repeats
    sout +=f"\n{indent}REPEAT {rep} " + "{"
    indent+="    "

    sout=do_round(sout)
    
    if p[0]>0:                        # depolarize just before measurement
        sout += f"\n{indent}DEPOLARIZE1({p[0]}) "
        for i0 in range(ndat,ntot):
            sout += f" {i0}"

    sout += f"\n{indent}MR"
    for i0 in range(ndat,ntot):
        sout += f" {i0}"
    
    # insert detector events comparing this round and last round ancilla measurements
    for i in range(nd):
        sout += f"\n{indent}DETECTOR"
        for j in range(len(dmatrix[i])):
            sout +=f" rec[{-na-1+dmatrix[i][j]}] rec[{-2*na-1+dmatrix[i][j]}]"

    indent="    "
    sout += f"\n{indent}"+"}"

    # depolarize each data qubit just before the final measurement
    if p[0]>0: 
        sout += f"\n{indent}DEPOLARIZE1({p[0]})"
        for i0 in range(ndat):
            sout += f" {i0}"

    sout += f"\n{indent}M"
    for i0 in range(ndat):
        sout += f" {i0}"

    # insert last round detector events
    #print(set(dlastmatrix[len(dlastmatrix)-1]))
    for i in range(hnd,nd):
    #for i in range(len(dlastmatrix)):
        sout += f"\n{indent}DETECTOR"
        for j in range(len(dlastmatrix[i])):
            sout += f" rec[{-ndat+dlastmatrix[i][j]}]" # data qubits of the stabilizer 
        #sout +=f" rec[{-ndat-na+i}] "         # ancilla qubit of stabilizers on 

    sout +=f"\n{indent}TICK"

    #logical operator
    for i in range(nl):
        sout += f"\n{indent}OBSERVABLE_INCLUDE({i})"
        for j in range(len(savelnums[i])):
            sout +=f" rec[{-ndat+int(savelnums[i][j])}]"

    sout += f"\n# end of the circuit\n"
    return sout

################################################################
