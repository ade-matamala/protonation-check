from biopandas.pdb import PandasPdb
import pandas as pd
import nglview as ngv
from collections import Counter 

model0 = PandasPdb().read_pdb('model0_chainID.pdb') #read the pdb file
#ngv.show_file('model0_chainID.pdb') #vizualization of the structure 
df = model0.df['ATOM'] 

LSU = ['A','C','E','G','I','K','M','O'] #large subunit 
SSU = ['B','D','F','H','J','L','N','P'] #small subunit
chains = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P'] 
PASR_name=['GLU','ASP','KCX','HID','HIE','HIP','ARG','LYS'] 

#protonation states of each PASR_name
acid = {'GLU-':['N','H','CA','HA','CB','HB2','HB3','CG','HG2','HG3','CD','OE1','OE2','C','O'], 
        'GLU0':['N','H','CA','HA','CB','HB2','HB3','CG','HG2','HG3','CD','OE1','OE2','HE2','C','O'],
        'KCX-':['N','CA','CB','CG','CD','CE','NZ','C','O','CX','OQ1','OQ2','H13','H14','H15','H16','H3','H5','H6','H20','H21','H22', 'H23'],
        'KCX0':['N','CA','CB','CG','CD','CE','NZ','C','O','CX','OQ1','OQ2','H13','H14','H15','H16','H3','H5','H6','H20','H21','H22','H23','H27'],
        'ASP-':['N','H','CA','HA','CB','HB2','HB3','CG','OD1','OD2','C','O'],
        'ASP0':['N','H','CA','HA','CB','HB2','HB3','CG','OD1','OD2','HD2','C','O']} 
basic = {'LYS+':['N','H','CA','HA','CB','HB2','HB3','CG','HG2','HG3','CD','HD2','HD3','CE','HE2','HE3','NZ','HZ1','HZ2','HZ3','C','O'], 
         'LYS0':['N','H','CA','HA','CB','HB2','HB3','CG','HG2','HG3','CD','HD2','HD3','CE','HE2','HE3','NZ','HZ1','HZ2','C','O'],
         'HIP+':['N','H','CA','HA','CB','HB2','HB3','CG','ND1','HD1','CE1','HE1','NE2','HE2','CD2','HD2','C','O'],
         'HIE0':['N','H','CA','HA','CB','HB2','HB3','CG','ND1','CE1','HE1','NE2','HE2','CD2','HD2','C','O'],
         'HID0':['N','H','CA','HA','CB','HB2','HB3','CG','ND1','HD1','CE1','HE1','NE2','CD2','HD2','C','O'],         
         'ARG+':['N','H','CA','HA','CB','HB2','HB3','CG','HG2','HG3','CD','HD2','HD3','NE','HE','CZ','NH1','1HH1','2HH1','NH2','1HH2','2HH2','C','O'],
         'ARG0':['N','H','CA','HA','CB','HB2','HB3','CG','HG2','HG3','CD','HD2','HD3','NE','HE','CZ','NH1','2HH1','NH2','1HH2','2HH2','C','O']}

acid_base = [] #A, B
charge = [] #1, 0, -1
name = [] 
id_number = [] 
chain_id = [] 

for chain in chains: 
    if chain in LSU: 
        for res in range(1,451):
            residue_df = df.loc[(df['chain_id']==chain)&(df['residue_number']==res)]
            residue_name = residue_df['residue_name'].unique()
            cntr = Counter(residue_df['atom_name'])
    
            if residue_name in ['GLU','ASP','KCX']: #acid
                if cntr == Counter(acid[residue_name[0]+'0']):
                    #print(chain,res,residue_name,'acid', 'neutral')
                    acid_base.append('A')
                    charge.append(0)
                    name.append(residue_name)
                    id_number.append(res)
                    chain_id.append(chain) 
                    
                    
                elif cntr == Counter(acid[residue_name[0]+'-']):
                    #print(chain,res,residue_name,'acid', 'negative')
                    acid_base.append('A')
                    charge.append(-1)
                    name.append(residue_name)
                    id_number.append(res)
                    chain_id.append(chain)
                    
                else:
                    print('Error_Acid')
                    
            if residue_name in ['LYS', 'ARG']: #basic
                if cntr == Counter(basic[residue_name[0]+'0']):
                    #print(chain,res,residue_name,'basic', 'neutral')
                    acid_base.append('B')
                    charge.append(0)
                    name.append(residue_name)
                    id_number.append(res)
                    chain_id.append(chain)
                    
                elif cntr == Counter(basic[residue_name[0]+'+']):
                    #print(chain,res,residue_name,'basic', 'positive')
                    acid_base.append('B')
                    charge.append(1)
                    name.append(residue_name)
                    id_number.append(res)
                    chain_id.append(chain)
                    
                else:
                    print('Error_Base')
                    
            if residue_name in ['HIE', 'HID', 'HIP']: #histidines
                if cntr == Counter(basic[residue_name[0]+'0']):
                    #print(chain,res,residue_name,'basic', 'neutral')
                    acid_base.append('B')
                    charge.append(0)
                    name.append(residue_name)
                    id_number.append(res)
                    chain_id.append(chain)

                    
                elif cntr == Counter(basic['HIP'+'+']):
                    #print(chain,res,residue_name,'basic', 'positive')
                    acid_base.append('B')
                    charge.append(1)
                    name.append(residue_name)
                    id_number.append(res)
                    chain_id.append(chain)
                    
                else:
                    print('Error_Base')
   
    if chain in SSU:
        for res in range(1,100):
            residue_df = df.loc[(df['chain_id']==chain)&(df['residue_number']==res)]
            residue_name = residue_df['residue_name'].unique()
            cntr = Counter(residue_df['atom_name'])
          
            if residue_name in ['GLU','ASP','KCX']: #acid
                if cntr == Counter(acid[residue_name[0]+'0']):
                    #print(chain,res,residue_name,'acid', 'neutral')
                    acid_base.append('A')
                    charge.append(0)
                    name.append(residue_name)
                    id_number.append(res)
                    chain_id.append(chain)
                    
                elif cntr == Counter(acid[residue_name[0]+'-']):
                    #print(chain,res,residue_name,'acid', 'negative')
                    acid_base.append('A')
                    charge.append(-1)
                    name.append(residue_name)
                    id_number.append(res)
                    chain_id.append(chain)
                    
                else:
                    print('Error_Acid')
                    
            if residue_name in ['LYS', 'ARG']: #basic
                if cntr == Counter(basic[residue_name[0]+'0']):
                    #print(chain,res,residue_name,'basic', 'neutral')
                    acid_base.append('B')
                    charge.append(0)
                    name.append(residue_name)
                    id_number.append(res)
                    chain_id.append(chain)

                    
                elif cntr == Counter(basic[residue_name[0]+'+']):
                    #print(chain,res,residue_name,'basic', 'positive')
                    acid_base.append('B')
                    charge.append(1)
                    name.append(residue_name)
                    id_number.append(res)
                    chain_id.append(chain)
                    
                else:
                    print('Error_Base')
                    
            if residue_name in ['HIE', 'HID', 'HIP']: #histidines
                if cntr == Counter(basic[residue_name[0]+'0']):
                    #print(chain,res,residue_name,'basic', 'neutral')
                    acid_base.append('B')
                    charge.append(0)
                    name.append(residue_name)
                    id_number.append(res)
                    chain_id.append(chain)
                    
                elif cntr == Counter(basic['HIP'+'+']):
                    #print(chain,res,residue_name,'basic', 'positive')
                    acid_base.append('B')
                    charge.append(1)
                    name.append(residue_name)
                    id_number.append(res)
                    chain_id.append(chain)
                    
                else:
                    print('Error_HIS')

d = {'chain_id': chain_id, 'res_id': id_number, 'res_name':name,'acid/base': acid_base, 'charge': charge} 
data = pd.DataFrame(data=d) #create a dataframe with the data
data.to_excel('protonation.xlsx') #save the data in an excel file
