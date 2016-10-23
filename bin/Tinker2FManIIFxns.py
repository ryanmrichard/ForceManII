import os
import sys
import math
import subprocess as sub
usage="Usage: Tinker2FManII.py <path/2/xyz> <path/2/Params> [<path/2/analyze>]"
ang2au=1.889725989
kcalmol2au=1/627.5096
deg2rad=math.pi/180
k2au=kcalmol2au/(ang2au**2)
anglek2au=kcalmol2au/(deg2rad**2)
bond_type="FManII::IntCoord_t::Bond"
angle_type="FManII::IntCoord_t::Angle"
k_param="FManII::Param_t::K"
r0_param="FManII::Param_t::r0"
intcoords={bond_type:[k_param,r0_param],
           angle_type:[k_param,r0_param]
}

def check(cond,msg):
    if not cond: raise RuntimeError(msg)

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
    mol_name=os.path.splitext(os.path.basename(xyz_file))[0]
    return xyz_file,ff_file,corr_answers,mol_name

def read_sys(xyz_file):
    carts=[];connect=[];param_num=[]
    with open(xyz_file,"r") as f:
        next(f)
        for line in f:
            tokenized=line.split()
            carts.append([float(i) for i in tokenized[2:5]])
            param_num.append(int(tokenized[5]))
            connect.append([int(i)-1 for i in tokenized[6:]])#Tinker starts at 1
    return carts,connect,param_num

def read_ff(ff_file):
    atom2tink={};parms={ct:{p:{} for p in pt} for ct,pt in intcoords.items()}
    with open(ff_file,"r") as f:
        for line in f:
            da_line=line.split()
            if len(da_line)<1:continue
            if da_line[0]=="atom":atom2tink[int(da_line[1])]=int(da_line[2])
            if da_line[0]=="bond":
                i=int(da_line[1]);j=int(da_line[2])
                pair=(i,j) if i<=j else (j,i)
                #Tinker bakes the 1/2 into k already
                parms[bond_type][k_param][pair]=2*float(da_line[3])*k2au
                parms[bond_type][r0_param][pair]=float(da_line[4])*ang2au
            if da_line[0]=="angle":
                i=int(da_line[1]);j=int(da_line[2]);k=int(da_line[3])
                triple=(i,j,k) if i<=j else (k,j,i)
                parms[angle_type][k_param][triple]=2*float(da_line[4])*anglek2au
                parms[angle_type][r0_param][triple]=float(da_line[5])*deg2rad
    return atom2tink,parms
            
def compute_bonds(carts,atom2tink,connect,param_num,parms):
    bonds=[];bond_k=[];bond_r0=[]
    for i,r in enumerate(carts):
        for j in connect[i]:
            if j<i: continue #Only count bonds once
            ip=atom2tink[param_num[i]];jp=atom2tink[param_num[j]]
            pair=(ip,jp) if ip<=jp else (jp,ip)
            if not pair in parms[bond_type][r0_param]:
                print("No parameters for bond: "+str(pair))
                continue 
            rij=math.sqrt(sum([(r[k]-carts[j][k])**2 for k in range(0,3)]))
            bonds.append(rij*ang2au-parms[bond_type][r0_param][pair])
            #bonds.append((rij)*ang2au)
            bond_k.append(parms[bond_type][k_param][pair])
            bond_r0.append(parms[bond_type][r0_param][pair])
    return bonds,bond_k,bond_r0

def compute_angles(carts,atom2tink,connect,param_num,parms):
    angles=[];angle_k=[];angle_r0=[]
    for i,ri in enumerate(carts):
        for j in connect[i]:
            rj=carts[j]
            for k in connect[j]:
                if k==i:continue
                ip=atom2tink[param_num[i]];jp=atom2tink[param_num[j]];
                kp=atom2tink[param_num[k]]
                triple=(ip,jp,kp) if ip<=kp else (kp,jp,ip)
                if triple not in parms[angle_type][k_param]:
                    print("No value for angle "+str((i,j,k))+" "+str(triple))
                    continue
                rk=carts[k]
                rij=[ri[q]-rj[q] for q in range(0,3)]
                rjk=[rj[q]-rk[q] for q in range(0,3)]
                magrij=math.sqrt(rij[0]**2+rij[1]**2+rij[2]**2)
                magrjk=math.sqrt(rjk[0]**2+rjk[1]**2+rjk[2]**2)
                costheta=(rij[0]*rjk[0]+rij[1]*rjk[1]+rij[2]*rjk[2])/(magrij*magrjk)
                print(costheta)
                angles.append(math.acos(costheta))
                angle_k.append(parms[angle_type][k_param][triple])
                angle_r0.append(parms[angle_type][r0_param][triple])
    return angles,angle_k,angle_r0

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

def print_types(f,mol_name,param_num,atom2tink):
    f.write("const FManII::AtomTypes "+mol_name.lower()+"_FF_types={\n")
    for p in param_num:f.write("{"+str(atom2tink[p])+",0},\n")
    f.write("};\n\n")

def print_carts(f,mol_name,carts):
    f.write("const std::vector<double> "+mol_name.lower()+"={\n")
    for i in carts[:-1]:
        f.write(str(i[0]*ang2au)+","+str(i[1]*ang2au)+","+str(i[2]*ang2au)+",\n")
    f.write(str(carts[-1][0]*ang2au)+","+str(carts[-1][1]*ang2au)+
            ","+str(carts[-1][2]*ang2au)+"};\n\n")

def print_parms(f,mol_name,parms):
    """Prints parameter map in initializer list format
    Theactual map is coordtype by parameter type by atom types"""
    var_name=mol_name.lower()+"_FF_params"
    f.write("const FManII::ParamTypes "+var_name+"={\n")#outer most map
    for ct,rest in parms.items():
        f.write("{"+ct+",{\n")#Start of pair 1st map and start of middle map
        for pt,rest1 in rest.items():
            f.write("   {"+pt+",{\n")#Start of pair middle map and start of last map
            for atoms,value in rest1.items():
                f.write("      {{")#Start of pair last map and start of atom vector
                for ai in atoms:f.write(str(ai)+",")
                f.write("},"+str(value)+"},\n")#end atom vector, end pair last map
            f.write("      }},\n")#end last map end pair middle map
        f.write("   }},\n")#end middle map end 1st pair
    f.write("};\n")#end outer map
