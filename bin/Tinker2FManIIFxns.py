import os
import sys
import math
import subprocess as sub
#General notes about Tinker
#bond r0 terms are in Angstrom
#bond k are in kcal/(mol Angstrom^2)
#angle r0 are in degrees
#angle k is in kcal/(mol radian^2)
#torsion amplitude per path in kcal/mol
#torsion phase shift in degree

usage="Usage: Tinker2FManII.py <path/2/xyz> <path/2/Params> [<path/2/analyze>]"
ang2au=1.889725989
kcalmol2au=1/627.5096
deg2rad=math.pi/180
k2au=kcalmol2au/(ang2au**2)
anglek2au=kcalmol2au
bond_type="FManII::IntCoord_t::BOND"
angle_type="FManII::IntCoord_t::ANGLE"
torsion_type="FManII::IntCoord_t::TORSION"
imp_type="FManII::IntCoord_t::IMPTORSION"
vdw_type="FManII::IntCoord_t::LENNARD_JONES"
charge_type="FManII::IntCoord_t::ELECTROSTATICS"
k_param="FManII::Param_t::K"
r0_param="FManII::Param_t::r0"
v_param="FManII::Param_t::amp"
phi_param="FManII::Param_t::phi"
n_param="FManII::Param_t::n"
v2_param="FManII::Param_t::amp2"
phi2_param="FManII::Param_t::phi2"
n2_param="FManII::Param_t::n2"
v3_param="FManII::Param_t::amp3"
phi3_param="FManII::Param_t::phi3"
n3_param="FManII::Param_t::n3"
q_param="FManII::Param_t::q"
sigma_param="FManII::Param_t::sigma"
epsilon_param="FManII::Param_t::epsilon"
rad_rule=""
rad_type=""
epsilon_rule=""
vdw_scale=""
chg_scale=""
epsilon_0=""
dielectric=""
intcoords={bond_type:[k_param,r0_param],
           angle_type:[k_param,r0_param],
           torsion_type:[v_param,phi_param,n_param,
                         v2_param,phi2_param,n2_param,
                         v3_param,phi3_param,n3_param
                         ],
           imp_type:[v_param,phi_param,n_param],
           vdw_type:[sigma_param,epsilon_param],
           charge_type:[q_param]
}

def check(cond,msg):
    if not cond: raise RuntimeError(msg)

def cross(v1,v2):
    return [v1[1]*v2[2]-v1[2]*v2[1],
            v1[2]*v2[0]-v1[0]*v2[2],
            v1[0]*v2[1]-v1[1]*v2[0]]

def dot(v1,v2):
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]

def setup():
    corr_answers={}
    check(len(sys.argv)>=3,usage)
    xyz_file=sys.argv[1];ff_file=sys.argv[2]
    check(os.path.isfile(xyz_file),"XYZ file does not exist")
    check(os.path.isfile(ff_file),"Param file does not exit")
    if len(sys.argv)==4:
        xyz_fp=os.path.abspath(sys.argv[1])
        ff_fp=os.path.abspath(sys.argv[2])
        anal_fp=os.path.abspath(sys.argv[3])
        check(os.path.isfile(anal_fp) and os.access(anal_fp,os.X_OK),
          "analyze does note exist or is not executable")
        p=sub.Popen([anal_fp,xyz_fp,ff_fp,"E"],stdout=sub.PIPE)
        for line in p.stdout:
            da_line=line.split()
            if len(da_line)<=1:continue
            if da_line[0].decode('UTF-8')=="Bond":
                corr_answers["bond"]=float(da_line[2].decode('UTF-8'))*kcalmol2au
            elif da_line[0].decode("UTF-8")=="Angle":
                corr_answers["angle"]=float(da_line[2].decode('UTF-8'))*kcalmol2au
            elif da_line[0].decode("UTF-8")=="Torsional":
                corr_answers["torsion"]=float(da_line[2].decode('UTF-8'))*kcalmol2au
            elif da_line[0].decode("UTF-8")=="Improper":
                corr_answers["improper"]=float(da_line[2].decode('UTF-8'))*kcalmol2au
            elif da_line[0].decode("UTF-8")=="Van":
                corr_answers["vdw"]=float(da_line[3].decode('UTF-8'))*kcalmol2au
            elif da_line[0].decode("UTF-8")=="Charge-Charge":
                corr_answers["electrostatics"]=float(da_line[1].decode('UTF-8'))*kcalmol2au
                
    mol_name=os.path.splitext(os.path.basename(xyz_file))[0]
    return xyz_file,ff_file,corr_answers,mol_name

def read_sys(xyz_file):
    carts=[];connect=[];param_num=[]
    with open(xyz_file,"r") as f:
        next(f)
        for line in f:
            tokenized=line.split()
            carts.append([float(i)*ang2au for i in tokenized[2:5]])
            param_num.append(int(tokenized[5]))
            connect.append([int(i)-1 for i in tokenized[6:]])#Tinker starts at 1
    return carts,connect,param_num

def read_ff(ff_file):
    atom2tink={};parms={ct:{p:{} for p in pt} for ct,pt in intcoords.items()}
    with open(ff_file,"r") as f:
        for line in f:
            da_line=line.split()
            if len(da_line)<1:continue
            if da_line[0]=="atom":
                atom2tink[int(da_line[1])]=int(da_line[2])
            if da_line[0]=="bond":
                i=int(da_line[1]);j=int(da_line[2])
                pair=(i,j) if i<=j else (j,i)
                #Tinker bakes the 1/2 into k already
                parms[bond_type][k_param][pair]=2*float(da_line[3])*k2au
                parms[bond_type][r0_param][pair]=float(da_line[4])*ang2au
            if da_line[0]=="angle":
                i=int(da_line[1]);j=int(da_line[2]);k=int(da_line[3])
                triple=(i,j,k) if i<=k else (k,j,i)
                parms[angle_type][k_param][triple]=2*float(da_line[4])*anglek2au
                parms[angle_type][r0_param][triple]=float(da_line[5])*deg2rad
            if da_line[0]=="torsion":
                i,j,k,l=int(da_line[1]),int(da_line[2]),int(da_line[3]),int(da_line[4])
                quad=(i,j,k,l) if (j<k or (j==k and i<l)) else (l,k,j,i)

                parms[torsion_type][v_param][quad]=float(da_line[5])*kcalmol2au
                parms[torsion_type][phi_param][quad]=float(da_line[6])*deg2rad
                parms[torsion_type][n_param][quad]=float(da_line[7])
                if len(da_line)>8:
                    parms[torsion_type][v2_param][quad]=float(da_line[8])*kcalmol2au
                    parms[torsion_type][phi2_param][quad]=float(da_line[9])*deg2rad
                    parms[torsion_type][n2_param][quad]=float(da_line[10])
                    if len(da_line)==14:
                        parms[torsion_type][v3_param][quad]=float(da_line[11])*kcalmol2au
                        parms[torsion_type][phi3_param][quad]=float(da_line[12])*deg2rad
                        parms[torsion_type][n3_param][quad]=float(da_line[13])
            if da_line[0]=="imptors":
                i,j,k,l=int(da_line[1]),int(da_line[2]),int(da_line[3]),int(da_line[4])
                jkp=[i,j,l]
                jkp.sort()
                quad=(jkp[0],k,jkp[1],jkp[2])
                parms[imp_type][v_param][quad]=float(da_line[5])*kcalmol2au
                parms[imp_type][phi_param][quad]=float(da_line[6])*deg2rad
                parms[imp_type][n_param][quad]=float(da_line[7])
            if da_line[0]=="charge":
                i=int(da_line[1])
                elem=(i,)
                parms[charge_type][q_param][elem]=float(da_line[2])
            if da_line[0]=="vdw":
                i=int(da_line[1])
                elem=(i,)
                parms[vdw_type][sigma_param][elem]=\
                float(da_line[2])*ang2au*(2.0 if rad_type=="R-MIN" else 1.0)
                parms[vdw_type][epsilon_param][elem]=float(da_line[3])*kcalmol2au
            if da_line[0]=="radiusrule":
                global rad_rule
                rad_rule=da_line[1]
            if da_line[0]=="radiustype":
                global rad_type
                rad_type=da_line[1]
            if da_line[0]=="epsilonrule":
                global epsilon_rule
                epsilon_rule=da_line[1]
            if da_line[0]=="vdw-14-scale":
                global vdw_scale
                vdw_scale=da_line[1]
            if da_line[0]=="chg-14-scale":
                global chg_scale
                chg_scale=da_line[1]
            if da_line[0]=="electric":
                global epsilon_0
                epsilon_0=da_line[1]
            if da_line[0]=="dielectric":
                global dielectric
                dielectric=da_line[1]
                
    return atom2tink,parms,rad_rule,rad_type,epsilon_rule,vdw_scale,chg_scale
            
def compute_bonds(carts,atom2tink,connect,param_num,parms):
    bonds=[];bond_k=[];bond_r0=[];
    for i,r in enumerate(carts):
        for j in connect[i]:
            if j<i: continue #Only count bonds once
            ip=atom2tink[param_num[i]];jp=atom2tink[param_num[j]]
            pair=(ip,jp) if ip<=jp else (jp,ip)
            if not pair in parms[bond_type][r0_param]:
                print("No parameters for bond: "+str(pair))
                continue 
            rij=math.sqrt(sum([(r[k]-carts[j][k])**2 for k in range(0,3)]))
            bonds.append(rij-parms[bond_type][r0_param][pair])
            #bonds.append((rij)*ang2au)
            bond_k.append(parms[bond_type][k_param][pair])
            bond_r0.append(parms[bond_type][r0_param][pair])
    return bonds,bond_k,bond_r0
def compute_pairs(carts,atom2tink,connect,param_num,parms):
    pairs,chg,vdw=[],[],[]
    all_12_pairs,all_13_pairs,all_14_pairs=set(),set(),set()
    for i,r in enumerate(carts):
        qi=parms[charge_type][q_param][(param_num[i],)]
        si=parms[vdw_type][sigma_param][(atom2tink[param_num[i]],)]
        ei=parms[vdw_type][epsilon_param][(atom2tink[param_num[i]],)]
        for j in connect[i]:
            for k in connect[j]:
                if k==i:continue
                for l in connect[k]:
                    if l==j:continue
                    if l<i: continue
                    if l in connect[i]:continue
                    all_14_pairs.add((i,l))
                if k<i:continue
                all_13_pairs.add((i,k))
            if j<i:continue
            all_12_pairs.add((i,j))
    for i,ri in enumerate(carts):
        qi=parms[charge_type][q_param][(param_num[i],)]
        si=parms[vdw_type][sigma_param][(atom2tink[param_num[i]],)]
        ei=parms[vdw_type][epsilon_param][(atom2tink[param_num[i]],)]
        for j in range(i+1,len(carts)):
            is12=(i,j) in all_12_pairs
            is13=(i,j) in all_13_pairs
            is14=(i,j) in all_14_pairs
            if is12 or is13:continue
            rij=[ri[x]-carts[j][x] for x in range(3)]
            magrij=math.sqrt(dot(rij,rij))
            qj=parms[charge_type][q_param][(param_num[j],)]
            sj=parms[vdw_type][sigma_param][(atom2tink[param_num[j]],)]
            ej=parms[vdw_type][epsilon_param][(atom2tink[param_num[j]],)]
            chg.append(qi*qj*(1/float(chg_scale) if is14 else 1.0))
            pairs.append(magrij)
            vdw.append((math.sqrt(ei*ej)*(1/float(vdw_scale) if is14 else 1.0),(si+sj)))
            #print((i+1,j+1),qi*qj/magrij/kcalmol2au)

    return pairs,chg,vdw
                    

def compute_angles(carts,atom2tink,connect,param_num,parms):
    angles=[];angle_k=[];angle_r0=[]
    for i,ri in enumerate(carts):
        for j in connect[i]:
            rj=carts[j]
            for k in connect[j]:
                if k<=i:continue
                ip=atom2tink[param_num[i]];jp=atom2tink[param_num[j]];
                kp=atom2tink[param_num[k]]
                triple=(ip,jp,kp) if ip<=kp else (kp,jp,ip)
                if triple not in parms[angle_type][k_param]:
                    print("No value for angle "+str((i,j,k))+" "+str(triple))
                    continue
                rk=carts[k]
                rij=[ri[q]-rj[q] for q in range(3)]
                rjk=[rk[q]-rj[q] for q in range(3)]
                magrij=math.sqrt(dot(rij,rij))
                magrjk=math.sqrt(dot(rjk,rjk))
                costheta=dot(rij,rjk)/(magrij*magrjk)
                angles.append(math.acos(costheta)-parms[angle_type][r0_param][triple])
                angle_k.append(parms[angle_type][k_param][triple])
                angle_r0.append(parms[angle_type][r0_param][triple])
    return angles,angle_k,angle_r0
def tor_common(ri,rj,rk,rl):
    r21=[ri[q]-rj[q] for q in range(3)]
    r23=[rj[q]-rk[q] for q in range(3)]
    r34=[rl[q]-rk[q] for q in range(3)]
    n1=cross(r21,r23)
    n2=cross(r34,r23)
    magr23=math.sqrt(dot(r23,r23))
    magn1=math.sqrt(dot(n1,n1))
    magn2=math.sqrt(dot(n2,n2))
    n1dotn2=dot(n1,n2)
    n1crossn2=cross(n1,n2)
    n1cn2dotn3=dot(n1crossn2,r23)
    cosphi=n1dotn2/(magn1*magn2)
    sinphi=n1cn2dotn3/(magn1*magn2*magr23)
    phi=math.atan2(sinphi,cosphi)
    return phi

def compute_torsions(carts,atom2tink,connect,param_num,parms):
    torsions=[];torsion_v=[];
    for i,ri in enumerate(carts):
        ip=atom2tink[param_num[i]]
        for j in connect[i]:
            rj=carts[j]
            jp=atom2tink[param_num[j]]
            for k in connect[j]:
                if k==i:continue
                if k<j:continue
                rk=carts[k]
                kp=atom2tink[param_num[k]]
                for l in connect[k]:
                    if l==j:continue
                    rl=carts[l]
                    r14=[rl[x]-ri[x] for x in range(3)]
                    lp=atom2tink[param_num[l]]
                    quad=(ip,jp,kp,lp) if (jp<kp or (jp==kp and ip<lp)) else (lp,kp,jp,ip)
                    phi=tor_common(ri,rj,rk,rl)
                    if quad not in parms[torsion_type][n_param]:
                        print(str(quad)+"is not in the params. atoms: "+str((i,j,k,l)))
                        continue
                    n=parms[torsion_type][n_param][quad]
                    gamma=parms[torsion_type][phi_param][quad]
                    v=parms[torsion_type][v_param][quad]
                    torsions.append(n*phi-gamma)
                    torsion_v.append(v)
                    if quad in parms[torsion_type][n2_param]:
                        n=parms[torsion_type][n2_param][quad]
                        gamma=parms[torsion_type][phi2_param][quad]
                        v=parms[torsion_type][v2_param][quad]
                        torsions.append(n*phi-gamma)
                        torsion_v.append(v)
                    if quad in parms[torsion_type][n3_param]:
                        n=parms[torsion_type][n3_param][quad]
                        gamma=parms[torsion_type][phi3_param][quad]
                        v=parms[torsion_type][v3_param][quad]
                        torsions.append(n*phi-gamma)
                        torsion_v.append(v)
    return torsions,torsion_v
def compute_imp_torsions(carts,atom2tink,connect,param_num,parms):
    torsions=[];torsion_v=[]
    for i,ri in enumerate(carts):
        ip=atom2tink[param_num[i]]
        for j in connect[i]:
            if len(connect[j])!=3:continue
            rj=carts[j]
            jp=atom2tink[param_num[j]]
            for k in connect[j]:
                if k<=i:continue
                rk=carts[k]
                kp=atom2tink[param_num[k]]
                for l in connect[j]:
                    if l<=k:continue
                    rl=carts[l]
                    lp=atom2tink[param_num[l]]
                    jkp=[ip,kp,lp]

                    rji=[rj[x]-ri[x] for x in range(3)]
                    rjk=[rj[x]-rk[x] for x in range(3)]
                    rjl=[rj[x]-rl[x] for x in range(3)]
                    magji=math.sqrt(dot(rji,rji))
                    magjk=math.sqrt(dot(rjk,rjk))
                    magjl=math.sqrt(dot(rjl,rjl))
                    da_order=[]
                    if ip==kp:
                        angle1=math.acos(dot(rjl,rji)/(magji*magjl))
                        angle2=math.acos(dot(rjl,rjk)/(magjk*magjl))
                        deg=(0,1) if angle1<angle2 else (1,0)
                        da_order=[deg[0],deg[1],2] if lp>ip else [2,deg[0],deg[1]]
                    elif ip==lp:
                        angle1=math.acos(dot(rjk,rji)/(magjk*magji))
                        angle2=math.acos(dot(rjk,rjl)/(magjk*magjl))
                        deg=(0,2) if angle1<angle2 else (2,0)
                        da_order=[deg[0],deg[1],1] if kp>ip else [1,deg[0],deg[1]]
                    elif kp==lp:
                        angle1=math.acos(dot(rji,rjk)/(magji*magjk))
                        angle2=math.acos(dot(rji,rjl)/(magji*magjl))
                        deg=(1,2) if angle1<angle2 else (2,1)
                        da_order=[deg[0],deg[1],0] if ip>kp else [0,deg[0],deg[1]]
                    else:
                        da_order=sorted(range(len(jkp)), key=lambda k: jkp[k])
                    quad=(jkp[da_order[0]],jp,jkp[da_order[1]],jkp[da_order[2]])
                    if quad not in parms[imp_type][n_param]:
                        print(str(quad)+"is not in the params. atoms: "+str((i,j,k,l)))
                        continue
                    rs=[ri,rk,rl]
                    phi=tor_common(rs[da_order[0]],rj,rs[da_order[1]],rs[da_order[2]])

                    n=parms[imp_type][n_param][quad]
                    gamma=parms[imp_type][phi_param][quad]
                    v=parms[imp_type][v_param][quad]

                    torsions.append(n*phi-gamma)
                    torsion_v.append(v)
    return torsions,torsion_v
def print_conns(f,mol_name,connect):
    f.write("const FManII::ConnData "+mol_name.lower()+"_conns={\n")
    for i in connect[:-1]:
        f.write("{")#Start set of connected atoms
        for j in i[:-1]:f.write(str(j)+",")
        f.write(str(i[-1])+"},\n")#end set of connected atoms
    f.write("{")#Last atom set that doesn't get a comma
    for j in connect[-1][:-1]:f.write(str(j)+",")
    f.write(str(connect[-1][-1])+"}};\n\n")

def print_sys_parms(f,mol_name,desc,array):
    f.write("const std::vector<double> "+mol_name.lower()+"_"+desc+"{\n")
    for i in array:
        f.write(str(i)+",\n")
    f.write("};\n")
    
def print_coord(f,mol_name,bonds_in,desc):
    f.write("const std::vector<double> "+mol_name.lower()+"_"+desc+"={\n")
    for i in bonds_in[:-1]:f.write(str(i)+",\n")
    f.write(str(bonds_in[-1])+"};\n\n")

def print_carts(f,mol_name,carts):
    f.write("const std::vector<double> "+mol_name.lower()+"={\n")
    for i in carts[:-1]:
        f.write(str(i[0])+","+str(i[1])+","+str(i[2])+",\n")
    f.write(str(carts[-1][0])+","+str(carts[-1][1])+
            ","+str(carts[-1][2])+"};\n\n")

