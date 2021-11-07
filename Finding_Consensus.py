#!/usr/bin/env python
# coding: utf-8

# In[80]:


find=open('rosalind_cons(2).txt')
stringg=""
lstStr=[]
d={}

#А ЭТО НОРМ
for line in find:
        if line.startswith('>'):
            key=line.rstrip()
            b=''
        else:
            b+=line.rstrip()
        d[key]=b       
        
for val in d.values():
    lstStr.append(val)

def prof_matrix(a):
    matr=[]
    for i in range(4):
        matr.append([0]*len(a[0]))
    prof_dct={"A": matr[0][:], "C": matr[1][:], "G": matr[2][:],"T": matr[3][:]}
    for word in a:
        for count, c in enumerate(word):
            prof_dct[c][count]+=1
        
    return prof_dct

def consensus(prof_matrix):
    keys=["A","T","G","C"]
    result_string=[]
    for i in range(len(prof_matrix[keys[0]])):
        fixed_max_count=0
        fixed_max_letter=None
        for k in keys:
            if prof_matrix[k][i]>fixed_max_count:
                fixed_max_count=prof_matrix[k][i]
                fixed_max_letter=k
        result_string.append(fixed_max_letter)
    return "".join(result_string)
print(consensus(prof_matrix(lstStr)))


   
find.close()
find_write=open("jj.txt","w")
find_write.write(consensus(prof_matrix(lstStr))+"\n")
for key,val in prof_matrix(lstStr).items():
    myStr=""
    val=" ".join([str(x) for x in val])
    myStr+=key+": "+" "+val
    find_write.write(myStr+"\n")
    print(myStr)
find_write.close()


# In[79]:


#СПИЗЖЕННЫЙ СКРИПТ, ДЕЛАЕТ ТОЖЕ САМОЕ
from Bio import SeqIO                      
sequences = []                             
handle = open('rosalind_cons.txt', 'r')
for record in SeqIO.parse(handle, 'fasta'):
     sequence = []                         
     for nt in record.seq:                 
          sequence.extend(nt)              
     sequences.append(sequence)

handle.close() 
profile = [[0]*len(sequences)]*4
import numpy                                                  
profile = numpy.zeros((4, len(sequences[0])), dtype=numpy.int)
for i,line in enumerate(sequences):                           
     for j, nt in enumerate(line):                            
          if nt == 'A':                                       
               profile[0][j] += 1                             
          elif nt == 'C':                                     
               profile[1][j] += 1                             
          elif nt == 'G':                                     
               profile[2][j] += 1                             
          elif nt == 'T':                                     
               profile[3][j] += 1    
consensus = ''                                                  
for A,C,G,T in zip(profile[0],profile[1],profile[2],profile[3]):
     if A >= C and A >= G and A >= T:                           
          consensus += 'A'                                      
     elif C >= A and C >= G and C >= T:                         
          consensus += 'C'                                      
     elif G >= A and G >= C and G >= T:                         
          consensus += 'G'                                      
     elif T >= A and T >= C and T >= G:                         
          consensus += 'T'        
            
print(consensus)
print('A: ' + ' '.join(str(e) for e in profile[0]))
print('C: ' + ' '.join(str(e) for e in profile[1]))
print('G: ' + ' '.join(str(e) for e in profile[2]))
print('T: ' + ' '.join(str(e) for e in profile[3]))


# In[ ]:




