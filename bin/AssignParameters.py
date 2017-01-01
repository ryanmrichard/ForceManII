from ParseFF import *

def compute_bonds(carts,atom2tink,connect,param_num,parms):
    bond_k,bond_r0=[],[]
    for i,r in enumerate(carts):
        for j in connect[i]:
            if j<i: continue #Only count bonds once
            ip,jp=atom2tink[param_num[i]],atom2tink[param_num[j]]
            pair=(ip,jp) if ip<=jp else (jp,ip)
            ffterm=(models["ho"],intcoords["bond"])
            if not pair in parms[ffterm][params["r0"]]:
                print("No parameters for bond: "+str(pair))
                continue
            bond_k.append(parms[ffterm][params["k"]][pair][0])
            bond_r0.append(parms[ffterm][params["r0"]][pair][0])
    return bond_k,bond_r0

def compute_angles(carts,atom2tink,connect,param_num,parms):
    angle_k,angle_r0=[],[]
    for i,ri in enumerate(carts):
        ip=atom2tink[param_num[i]]
        for j in connect[i]:
            jp=atom2tink[param_num[j]]
            for k in connect[j]:
                if k<=i:continue
                kp=atom2tink[param_num[k]]
                triple=(ip,jp,kp) if ip<=kp else (kp,jp,ip)
                ffterm=(models["ho"],intcoords["angle"])
                if triple not in parms[ffterm][params["k"]]:
                    print("No value for angle "+str((i,j,k))+" "+str(triple))
                    continue
                angle_k.append(parms[ffterm][params["k"]][triple][0])
                angle_r0.append(parms[ffterm][params["r0"]][triple][0])
    return angle_k,angle_r0

def compute_torsions(carts,atom2tink,connect,param_num,parms):
    torsion_v,torsion_phi,torsion_n=[],[],[]
    for i,ri in enumerate(carts):
        ip=atom2tink[param_num[i]]
        for j in connect[i]:
            jp=atom2tink[param_num[j]]
            for k in connect[j]:
                if k==i:continue
                if k<j:continue
                kp=atom2tink[param_num[k]]
                for l in connect[k]:
                    if l==j:continue
                    lp=atom2tink[param_num[l]]
                    quad=(ip,jp,kp,lp) if (jp<kp or (jp==kp and ip<lp)) else (lp,kp,jp,ip)
                    ffterm=(models["fs"],intcoords["torsion"])
                    if quad not in parms[ffterm][params["n"]]:
                        print(str(quad)+"is not in the params. atoms: "+str((i,j,k,l)))
                        continue
                    vs=parms[ffterm][params["v"]][quad]
                    phis=parms[ffterm][params["phi"]][quad]
                    ns=parms[ffterm][params["n"]][quad]
                    for m in range(3):
                        torsion_v.append(0.0 if len(vs)<=m else vs[m])
                        torsion_phi.append(0.0 if len(phis)<=m else phis[m])
                        torsion_n.append(0.0 if len(ns)<=m else ns[m])
    return torsion_v,torsion_phi,torsion_n

def compute_imp_torsions(carts,atom2tink,connect,param_num,parms):
    torsion_v,torsion_phi,torsion_n=[],[],[]
    for i,ri in enumerate(carts):
        ip=atom2tink[param_num[i]]
        for j in connect[i]:
            if len(connect[j])!=3:continue
            jp=atom2tink[param_num[j]]
            for k in connect[j]:
                if k<=i:continue
                kp=atom2tink[param_num[k]]
                for l in connect[j]:
                    if l<=k:continue
                    lp=atom2tink[param_num[l]]
                    jkp=[ip,kp,lp]
                    da_order=sorted(range(len(jkp)), key=lambda k: jkp[k])
                    quad=(jkp[da_order[0]],jp,jkp[da_order[1]],jkp[da_order[2]])
                    ffterm=(models["fs"],intcoords["imp"])
                    if quad not in parms[ffterm][params["n"]]:
                        print(str(quad)+"is not in the params. atoms: "+str((i,j,k,l)))
                        continue
                    torsion_n.append(parms[ffterm][params["n"]][quad][0])
                    torsion_phi.append(parms[ffterm][params["phi"]][quad][0])
                    torsion_v.append(parms[ffterm][params["v"]][quad][0])
    return torsion_v,torsion_phi,torsion_n

def compute_pairs(carts,atom2tink,connect,param_num,parms):
    ps=[{"pair14":[],"pair":[]} for i in range(3)]
    all_12_pairs,all_13_pairs,all_14_pairs=set(),set(),set()
    for i,r in enumerate(carts):
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
        for j in range(i+1,len(carts)):
            is12=(i,j) in all_12_pairs
            is13=(i,j) in all_13_pairs
            is14=(i,j) in all_14_pairs
            if is12 or is13:continue
            itype="pair14" if is14 else "pair"
            ctype=(models["cl"],intcoords[itype])
            vtype=(models["lj"],intcoords[itype])
            qi=parms[ctype][params["q"]][(param_num[i],)][0]
            si=parms[vtype][params["sigma"]][(atom2tink[param_num[i]],)][0]
            ei=parms[vtype][params["epsilon"]][(atom2tink[param_num[i]],)][0]
            qj=parms[ctype][params["q"]][(param_num[j],)][0]
            sj=parms[vtype][params["sigma"]][(atom2tink[param_num[j]],)][0]
            ej=parms[vtype][params["epsilon"]][(atom2tink[param_num[j]],)][0]
            ps[0][itype].append(qi*qj)
            ps[1][itype].append(math.sqrt(ei*ej))
            ps[2][itype].append(0.5*(si+sj))
    return ps[0],ps[1],ps[2]


######################## The main script is below ##############################

def main():
    if len(sys.argv) !=3:
        raise RuntimeError("Usage: python3 AssignParameters.py <xyz_file> <param_file>")
    xyz_file,ff_file=sys.argv[1],sys.argv[2]
    if not os.path.isfile(ff_file):
        raise RuntimeError("Parameter file does not exist")
    if not os.path.isfile(xyz_file):
        raise RuntimeError("XYZ file does not exist")
    carts,connect,param_num=read_sys(xyz_file)
    ff=read_ff(ff_file)
    ff_name=os.path.splitext(ff_file)[0]
    mol_name=os.path.splitext(xyz_file)[0]
    foundparams={}
    for key,value in ff.params.items():
        foundparams[key]={}
    k,r0=compute_bonds(carts,ff.type2class,connect,param_num,ff.params)
    ffterm=(models["ho"],intcoords["bond"])
    foundparams[ffterm][params["k"]]=k
    foundparams[ffterm][params["r0"]]=r0

    ffterm=(models["ho"],intcoords["angle"])
    k,r0=compute_angles(carts,ff.type2class,connect,param_num,ff.params)
    foundparams[ffterm][params["k"]]=k
    foundparams[ffterm][params["r0"]]=r0

    ffterm=(models["fs"],intcoords["torsion"])
    v,phi,n=compute_torsions(carts,ff.type2class,connect,param_num,ff.params)
    foundparams[ffterm][params["v"]]=v
    foundparams[ffterm][params["phi"]]=phi
    foundparams[ffterm][params["n"]]=n

    ffterm=(models["fs"],intcoords["imp"])
    v,phi,n=compute_imp_torsions(carts,ff.type2class,connect,param_num,ff.params)
    foundparams[ffterm][params["v"]]=v
    foundparams[ffterm][params["phi"]]=phi
    foundparams[ffterm][params["n"]]=n


    qs,es,ss=compute_pairs(carts,ff.type2class,connect,param_num,ff.params)
    for itype in ["pair14","pair"]:
        ctype=(models["cl"],intcoords[itype])
        vtype=(models["lj"],intcoords[itype])
        foundparams[ctype][params["q"]]=qs[itype]
        foundparams[vtype][params["epsilon"]]=es[itype]
        foundparams[vtype][params["sigma"]]=ss[itype]

    f=open(mol_name+"_params.hpp","w")
    f.write("//This file is autogenerated from AssignParameters.py\n\n")
    f.write("//If you value your sanity do not try to manually edit it!!!\n\n")
    f.write("#include<ForceManII/FManII.hpp>\n")
    f.write("extern const FManII::ParamSet "+mol_name+"_params;\n")
    f.close()

    f=open(mol_name+"_params.cpp","w")
    f.write("//This file is autogenerated from AssignParameters.py\n\n")
    f.write("//If you value your sanity do not try to manually edit it!!!\n\n")
    f.write("#include<ForceManII/FManII.hpp>\n")
    f.write("#include \""+mol_name+"_params.hpp\"\n")
    f.write("const FManII::ParamSet "+mol_name+"_params={\n")#Open map 1
    for key,val in foundparams.items():
        f.write("{{"+key[0]+","+key[1]+"},{\n")#Open map 1 item and open map 2
        for key2,val2 in val.items():
            f.write("{"+key2+",{\n")#Open map 2 item and open vector
            for ti in val2:
                f.write(str(ti)+",\n")
            f.write("}},\n")#Close vector close map 2 item
        f.write("}},\n")#Close map 2 close map 1 item
    f.write("};\n")#Close map 1
    f.close()



if __name__ == "__main__":
    main()
