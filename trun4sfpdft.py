#!/usr/bin/env python3

import numpy,h5py,sys
import h5py


navg = False
file2 = False
rohf = False

if len(sys.argv) > 2:
    try:
        nroot = int(sys.argv[2])
    except:
        if sys.argv[2] == 'avg':
            navg = True
        elif 'rohf' in sys.argv[2]:
            rohf = True
            hname_rohf = sys.argv[2]
        elif '.h5' in sys.argv[2]:
            file2 = True
            hname2 = sys.argv[2]
            
        else:
            print('Specify ROOT')
            sys.exit()
else:
    nroot = 0

    
if len(sys.argv) > 3:
    try:
        nroot = int(sys.argv[3])
    except:
        if sys.argv[3] == 'avg':
            navg = True
        else:
            s = sys.argv[3].split('.')
            if len(s) > 1:
                fname = str(sys.argv[3])
            else:
                fname = str(sys.argv[3])+'.Orb'
else:
    fname = 'orbital.Orb'

if len(sys.argv) > 4:
    s = sys.argv[4].split('.')
    if len(s) > 1:
        fname = str(sys.argv[4])
    else:
        fname = str(sys.argv[4])+'.Orb'
else:
    fname = 'orbital.Orb'
        

hname = sys.argv[1]

file3 = False

# read h5

f1 = h5py.File(hname,'r')

nbas = f1.attrs['NBAS'][0]
nsym = f1.attrs['NSYM']
dens_ = numpy.asarray(f1['DENSITY_MATRIX'])
mocc = numpy.asarray(f1['MO_OCCUPATIONS'])

if not rohf:
    print(' Using MOs from file ',hname)
    mo_coeff =  numpy.asarray(f1['MO_VECTORS'])
else:
    f_rohf = h5py.File(hname_rohf)
    if not nbas == f1.attrs['NBAS'][0]:
        print(' NBAS not the same with ROHF hdf5 file')
        sys.exit()

    print(' Using MOs from file ',hname_rohf)
    mo_coeff = numpy.asarray(f_rohf['MO_VECTORS'])
mo_coeff = mo_coeff.reshape((nbas,nbas),order='F')
mo_energy = numpy.asarray(f1['MO_ENERGIES'])
mo_type = numpy.asarray(f1['MO_TYPEINDICES'])
mo_type = mo_type.astype('str')
frozen = numpy.count_nonzero(mo_type=='F')
nras1 = numpy.count_nonzero(mo_type=='1')
nras2 = numpy.count_nonzero(mo_type=='2')
nras3 = numpy.count_nonzero(mo_type=='3')
rdm_type = mo_type[mo_type != 'F']


MO2 = False
# file2
if file2:
    print(' Using densities from 2 different files')
    f2 = h5py.File(hname2,'r')
    dens2 = numpy.asarray(f2['DENSITY_MATRIX'])
    mo_coeff2 =  numpy.asarray(f2['MO_VECTORS'])
    mo_coeff2 = mo_coeff2.reshape((nbas,nbas),order='F')
    mo_type2 = numpy.asarray(f2['MO_TYPEINDICES'])
    mo_type2 = mo_type.astype('str')
    mocc2 = numpy.asarray(f2['MO_OCCUPATIONS'])
    mo_energy2 = numpy.asarray(f2['MO_ENERGIES'])

    if not numpy.array_equal(mo_type,mo_type2):
        print(' RAS configuration not equal in the two files!')
        sys.exit()
    
    MO2 = True        
    
    # dens1 + den2
    ndens = len(dens_) + len(dens2)
    adens = dens_[0] + dens2[0]
    adens /= 2.0
apple = False

if not apple:
    if navg:
        if not file2:
            ndens = len(dens_)
            adens = sum(dens_)/float(ndens)
    else:
        adens = dens_[nroot]
    rdm_ras1 = adens[:nras1,:nras1]
    rdm_ras2 = adens[nras1:nras1+nras2,
                     nras1:nras1+nras2]
    rdm_ras3 = adens[nras1+nras2:,
                     nras1+nras2:]
    
    mo_ras1 = mo_coeff[:,mo_type == '1']
    mo_ras2 = mo_coeff[:,mo_type == '2']
    mo_ras3 = mo_coeff[:,mo_type == '3']

    ne1, nv1 = numpy.linalg.eigh(rdm_ras1)
    idx1 = numpy.argsort(-ne1)
    nat1 = numpy.dot(mo_ras1, nv1)

    ne2, nv2 = numpy.linalg.eigh(rdm_ras2)
    idx2 = numpy.argsort(-ne2)
    nat2 = numpy.dot(mo_ras2, nv2)

    ne3, nv3 = numpy.linalg.eigh(rdm_ras3)
    idx3 = numpy.argsort(-ne3)
    nat3 = numpy.dot(mo_ras3, nv3)

    natorb = numpy.zeros((nbas,nbas))
    natorb[:,mo_type == 'F'] = mo_coeff[:,mo_type == 'F']
    natorb[:,mo_type == '1'] = nat1[:,idx1]
    natorb[:,mo_type == '2'] = nat2[:,idx2]
    natorb[:,mo_type == '3'] = nat3[:,idx3]
    if MO2:
       natorba = numpy.zeros((nbas,nbas))
       mo_ras1a = mo_coeff2[:,mo_type == '1']
       mo_ras2a = mo_coeff2[:,mo_type == '2']
       mo_ras3a = mo_coeff2[:,mo_type == '3']

       nat1a = numpy.dot(mo_ras1a, nv1)
       nat2a = numpy.dot(mo_ras2a, nv2)
       nat3a = numpy.dot(mo_ras3a, nv3)

       natorba[:, mo_type == 'F'] = mo_coeff2[:, mo_type == 'F']
       natorba[:, mo_type == '1'] = nat1a[:, idx1]
       natorba[:, mo_type == '2'] = nat2a[:, idx2]
       natorba[:, mo_type == '3'] = nat3a[:, idx3]

    #tmpmo = mo_coeff[:,mo_type != 'F']
    #tmpe, tmpv = numpy.linalg.eigh(adens)
    #tmpidx = numpy.argsort(-tmpe)
    #tmpnat = numpy.dot(tmpmo, tmpv)
    #natorb[:,mo_type != 'F'] = tmpnat[:,tmpidx]    


    # print nat occ
    print(' Natural orbital occupation:')
    print(' Ras1 :',end='')
    for i in ne1[idx1]:
        print('{:>7.4f} '.format(i),end='')
    print()
    print(' Ras2 :',end='')
    for i in ne2[idx2]:
        print('{:>7.4f} '.format(i),end='')
    print()
    print(' Ras3 :',end='')
    for i in ne3[idx3]:
        print('{:>7.4f} '.format(i),end='')
    print()
    print()
    print(' RAS1 orb with occ < 1.998  : ',sum(ne1 < 1.998))
    print(' RAS3 orb with occ > 0.002  : ',sum(ne3 > 0.002))
    print()
    

    print(' RAS1 orb with occ < 1.999  : ',sum(ne1 < 1.999))
    print(' RAS3 orb with occ > 0.001  : ',sum(ne3 > 0.001))
    print()
    print(' RAS1 orb with occ < 1.9992  : ',sum(ne1 < 1.9992))
    print(' RAS3 orb with occ > 0.0008  : ',sum(ne3 > 0.0008))
    print()
    
    print(' RAS1 orb with occ < 1.9994  : ',sum(ne1 < 1.9994))
    print(' RAS3 orb with occ > 0.0006  : ',sum(ne3 > 0.0006))
    print()
    
    print(' RAS1 orb with occ < 1.9996  : ',sum(ne1 < 1.9996))
    print(' RAS3 orb with occ > 0.0004  : ',sum(ne3 > 0.0004))
    print()
    
    #print(' RAS1 orb with occ < 1.998  : ',sum(ne1 < 1.998))
    #print(' RAS3 orb with occ > 0.002  : ',sum(ne3 > 0.002))
    #print()

    print(' RAS1 orb with occ < 1.9998 : ',sum(ne1 < 1.9998))
    print(' RAS3 orb with occ > 0.0002 : ',sum(ne3 > 0.0002))
    print()

    print(' RAS1 orb with occ < 1.9999 : ',sum(ne1 < 1.99990e0))
    print(' RAS3 orb with occ > 0.0001 : ',sum(ne3 > 0.00010e0))
    print()

    print(' RAS1 orb with occ < 1.99992 : ',sum(ne1 < 1.99992))
    print(' RAS3 orb with occ > 0.00008 : ',sum(ne3 > 0.00008))
    print()

    print(' RAS1 orb with occ < 1.99998 : ',sum(ne1 < 1.99998))
    print(' RAS3 orb with occ > 0.00002 : ',sum(ne3 > 0.00002))
    print()

    print(' RAS1 orb with occ < 1.99999 : ',sum(ne1 < 1.99999))
    print(' RAS3 orb with occ > 0.00001 : ',sum(ne3 > 0.00001))
    print()



# ----
# Print molcas

# ----
if MO2:
   hs = hname.split('.')
   fs = fname.split('.')
   fnamea = fs[0] + hs[-1] + '.' + fs[-1]
else:
   fnamea = fname
   
w = open(fnamea, 'w')
w.write('#INPORB 2.2\n')
w.write('#INFO\n')
w.write('*blabla\n')
w.write('0   '+str(nsym)+'  0\n')
w.write(str(nbas)+'\n')
w.write(str(nbas)+'\n')

# Print orbitals
w.write('#ORB\n')
for i,j in enumerate(mo_coeff):
    n = j.shape[0]
    w.write('* ORBITAL '+str(1)+' '+str(i+1)+'\n')
    cout = 0
    for x in range(0,n,5):
        cout += 5 if n - x > 5 else n - x
        w.write('')
        for ic_ in range(x,cout):
            w.write('{:>22.14e}'.format(natorb[:,i][ic_]).upper())
        w.write('\n')

# Print occupation        
w.write('#OCC\n')
w.write('* OCCUPATION NUMBERS\n')
n = mocc.shape[0]
cout = 0
for x in range(0,n,5):
    cout += 5 if n - x > 5 else n - x
    for ic_ in range(x,cout):
        w.write('{:>22.14e}'.format(mocc[ic_]).upper())
    w.write('\n')
w.write('#OCHR\n')
w.write('* OCCUPATION NUMBERS (HUMAN-READABLE)\n')
for ic_ in mocc:
    w.write(' {:>7.4f} '.format(ic_))
w.write('\n')

# Print 1e- energies
w.write('#ONE\n')
w.write('* ONE ELECTRON ENERGIES\n')
w.write(' ')
for ic_ in mo_energy:
    w.write('{:>8.4e} '.format(ic_).upper())
w.write('\n')

w.write('#INDEX\n')
w.write('* 1234567890\n')
w.write('0 ')
for i in mo_type:
    w.write(i.lower())
w.write('\n')
w.close()
# ---


if MO2:
    hs = hname2.split('.')
    fs = fname.split('.')
    fnameb = fs[0] + hs[-1] + '.' + fs[-1]

    w = open(fnameb, 'w')
    w.write('#INPORB 2.2\n')
    w.write('#INFO\n')
    w.write('*blabla\n')
    w.write('0   '+str(nsym)+'  0\n')
    w.write(str(nbas)+'\n')
    w.write(str(nbas)+'\n')
    
    # Print orbitals
    w.write('#ORB\n')
    for i,j in enumerate(mo_coeff2):
        n = j.shape[0]
        w.write('* ORBITAL '+str(1)+' '+str(i+1)+'\n')
        cout = 0
        for x in range(0,n,5):
            cout += 5 if n - x > 5 else n - x
            w.write('')
            for ic_ in range(x,cout):
                w.write('{:>22.14e}'.format(natorba[:,i][ic_]).upper())
            w.write('\n')
    
    # Print occupation        
    w.write('#OCC\n')
    w.write('* OCCUPATION NUMBERS\n')
    n = mocc.shape[0]
    cout = 0
    for x in range(0,n,5):
        cout += 5 if n - x > 5 else n - x
        for ic_ in range(x,cout):
            w.write('{:>22.14e}'.format(mocc2[ic_]).upper())
        w.write('\n')
    w.write('#OCHR\n')
    w.write('* OCCUPATION NUMBERS (HUMAN-READABLE)\n')
    for ic_ in mocc2:
        w.write(' {:>7.4f} '.format(ic_))
    w.write('\n')
    
    # Print 1e- energies
    w.write('#ONE\n')
    w.write('* ONE ELECTRON ENERGIES\n')
    w.write(' ')
    for ic_ in mo_energy2:
        w.write('{:>8.4e} '.format(ic_).upper())
    w.write('\n')
    
    w.write('#INDEX\n')
    w.write('* 1234567890\n')
    w.write('0 ')
    for i in mo_type2:
        w.write(i.lower())
    w.write('\n')
    w.close()    
    # ---
    # ----
                                           

