#!/usr/bin/env python
# coding: utf-8

# In[1]:


## loading libraries

import os
import os.path
from os import listdir
from os.path import isfile, join
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statistics


# In[3]:


##listing the name of sequencing files
##cwd is the folder where the sequence files (txt) are
##cwd = os.getcwd()+"/txtfiles/"
seq_files = [f for f in listdir(cwd) if isfile(join(cwd, f))]


# In[6]:


## starting the pipeline

## looping over sequence files
for filename in ['QG exo1.txt']:

## skipping non sequence files
    if filename == '.jpg.png' or filename == '.csv' or filename == '.DS_Store':
        continue
    else:
## reading sequence files and excluding the empty or comment lines and giving unique number (line number) to each deletion event 
        with open(cwd+filename) as isolate:
            seqsdata = isolate.readlines()
            i=1
            sdas = []
            incs = []
            incst= []
            seqsdatanum = []
            compss =[]
            for sda in seqsdata:
                sdat= str(i)+sda
                seqsdatanum.append(sdat)
                incs.append(seqsdatanum.index(sdat)+1)
                i +=1
                if "#" in sdat or len(sdat) <50:
                    seqsdatanum.index(sdat)
                else:
                    incst.append(seqsdatanum.index(sdat)+1)
            for sda in seqsdata:
                if "#" in sda or len(sda) == 1:
                    sdas.append(sda)
            for ins in sdas:
                seqsdata.remove(ins)

## the first line in each txt file should be the reference (no deletion) sequence
        refseq = seqsdata[0]

## looping over each line (isolate) to identify the start and end sites of deletion as well as coordinates of flanking regions
        seqinfo = {}
        i = 0
        for seqsdataitem in seqsdata:
            f = seqsdataitem.find('-')
            rf = seqsdataitem.rfind('-')
            if f == -1:
                f = 0
                rf = len(seqsdata[0])
            bef = seqsdataitem[f-10:f].upper()
            aft = seqsdataitem[rf+1:rf+11].upper()
            rightbef = seqsdata[0][rf-9:rf+1].upper()
            leftafter = seqsdata[0][f:f+10].upper()
            refbef = seqsdata[0][f-10:f].upper()
            refaft = seqsdata[0][rf+1:rf+11].upper()
            refbefdelS = seqsdata[0][f+1:f+9].upper()
            refaftdelE = seqsdata[0][rf-6:rf+2].upper()
            
            comp1 = seqsdataitem[1:f].upper()
            comp2= seqsdataitem[rf+1:rf+20].upper()
            comp = str(comp1)+str(comp2)
            compss.append(comp)
            
## calculating length of deletion
            if i ==0:
                Dellength= 0
            else:
                Dellength = rf+1-f
                
            seqinfo[str(i)] = {'seqnum':i, 'start': f+1, 'end':rf+1, 'ref before':refbef, 'before':bef, 'ref leftafter':leftafter, 'ref after':refaft, 'after': aft, 'ref rightbefore':rightbef,'ref Start':refbefdelS, 'ref End':refaftdelE, 'Dellength':Dellength}
            i = i+1

## creating dataframe to contain retrieved info
        seqdf = pd.DataFrame.from_dict(seqinfo).transpose()[['seqnum', 'ref before', 'before','ref leftafter','start', 'Dellength','end', 'ref rightbefore','after', 'ref after']]

## ids corresponds to the line number in the sequence file
        ids = incst
        comps = compss
    
## searching for homology within 10 bp upstream and 10 bp downstream of the deletion junction
        matstrdicLARA={}
    
        for num in seqdf.iloc[1:,0]:
            LA = seqdf.iloc[num,3]
            RA= seqdf.iloc[num,8]
            matstr = ""
            for n in range(len(LA)):
                if LA[n] == RA[n]:
                    matstr += LA[n]
                else:
                    matstr += "-"
            matstrdicLARA['0'] = ""
            matstrdicLARA[num] = matstr

        matstrdicbefref={}

        for num in seqdf.iloc[1:,0]:
            bef = seqdf.iloc[num,2]
            refbef= seqdf.iloc[num,1]
            matstr = ""
            for n in range(len(bef)):
                if bef[n] == refbef[n]:
                    matstr += bef[n]
                else:
                    matstr += "-"
            matstrdicbefref['0'] = ""
            matstrdicbefref[num] = matstr

        matstrdicaftref={}

        for num in seqdf.iloc[1:,0]:
            aft = seqdf.iloc[num,8]
            refaft= seqdf.iloc[num,9]
            matstr = ""
            for n in range(len(aft)):
                if aft[n] == refaft[n]:
                    matstr += aft[n]
                else:
                    matstr += "-"
            matstrdicaftref['0'] = ""
            matstrdicaftref[num] = matstr

        matstrdicLBRB={}

        for num in seqdf.iloc[1:,0]:
            LB = seqdf.iloc[num,2]
            RB= seqdf.iloc[num,7]
            matstr = ""
            for n in range(len(bef)):
                if LB[n] == RB[n]:
                    matstr += LB[n]
                else:
                    matstr += "-"
            matstrdicLBRB['0'] = ""
            matstrdicLBRB[num] = matstr

## organizing the holomogly search results in a data frame
        seqdmatch = {'up ref vs sample':matstrdicaftref, 'down after vs up after':matstrdicLARA, 'down before vs up before':matstrdicLBRB, 'down ref vs sample':matstrdicbefref }
        seqdfmat = pd.DataFrame.from_dict(seqdmatch)
        
## adding the holomogly search results to the previously created data frame
        seqdf['up ref vs sample'] = list(seqdfmat['up ref vs sample'])
        seqdf['down after vs up after'] = list(seqdfmat['down after vs up after'])
        seqdf['down before vs up before'] = list(seqdfmat['down before vs up before'])
        seqdf['down ref vs sample'] = list(seqdfmat['down ref vs sample'])
        seqdf["counts"] = 1
        seqdf["staend"] = 1
        lissta = list(seqdf.iloc[:, 4])
        lisend = list(seqdf.iloc[:, 6])
        for ind in seqdf.index:
            sta= seqdf.loc[seqdf.index== ind, 'start']
            end = seqdf.loc[seqdf.index== ind, 'end']
            seqdf.loc[seqdf.index== ind, 'staend'] = str(int(sta))+str(int(end))
 
        seqdf.reset_index(drop=True,inplace=True)
        
## identifying the perfect microhology flanking deletion junctions
        ststs =[]
        for stst in seqdf['down after vs up after']:
            ch = stst.split('-')[0]
            if len(ch) >0:
                ststs.append(ch)
            else:
                ststs.append("")

        seqdf['down after vs up after homology'] = ststs
        
        ststss =[]
        for stst in seqdf['down before vs up before']:
            ch = stst.split('-')[-1]
            if len(ch) >0:
                ststss.append(ch)
            else:
                ststss.append("")
                
        seqdf['down before vs up before homology'] = ststss
        
        
        ststcom = []
        for ststnu in range(len(ststs)):
            ststcom.append(ststs[ststnu] + ststss[ststnu])
        
        seqdf['homology'] = ststcom

## adding the identified microhomology sequence to the data frame      
        for s in seqdf['down before vs up before homology']:
            if len(s) >1:
                if s[0] == '-':
                    seqdf['down before vs up before homology'][list(seqdf['down before vs up before homology']).index(s)]= seqdf['down before vs up before homology'][list(seqdf['down before vs up before homology']).index(s)].replace("-", " ", 1 )

## indentifying repeated deletion events
        seqdf['id'] = list(ids)
        seqdf['comp'] = list(comps)
        

        seqdf.index = seqdf['id']
        seqdf['repp'] = seqdf['id']
        delsam = {}        
        indexx = {}
        for dell in seqdf['Dellength'][seqdf['Dellength'].duplicated()]:
            delfil= seqdf['Dellength']==dell
            iddd= seqdf['id'][delfil]
            if len(iddd) > 1:
                startgroup = seqdf['start'][delfil]
                for startitem in startgroup:
                    delsam[str(list(seqdf['id'][delfil][abs(startgroup - startitem) < 15].values))]= seqdf['comp'][delfil][abs(startgroup - startitem) < 15] 
        i =0
        for i in delsam.keys():
            dd = delsam[i]
            ii =0
            indexxs={}
            for ii in dd.index:
                iii=0
                lenn=[]
                for iii in dd.index:
                    lenn.append(len(dd.loc[iii]))
                sml= min(lenn)
                dds = str(dd.loc[ii])[0:sml]
                indexxs[ii]= dds
            indexx[i]= indexxs
  
        filllslls = {}
        for gritkey in [*indexx]:
            i =0
            seit =[]
            for key, value in indexx[gritkey].items():
                seit.append(value)
            itmmseries = pd.Series([*indexx[gritkey]]) 
            seitlst = []
            for seitt in seit:
                seitlst.append(pd.Series(seit) == seitt)
            for st in seitlst:
                if len(itmmseries[st].tolist()) < 2:
                    continue
                else:
                    filllslls[str(gritkey)+str(i)] = itmmseries[st]
                    i =i+1
        
        groups =[]
        nuac = []
        for key, value in filllslls.items():
            if list(value) not in groups:
                groups.append(list(value))
                nuac.append(len(list(value)))
        i = 0
        seqdf['Occurrence'] = np.repeat(1, len(seqdf['id']))
        for groupitem in groups:
            for groupitemitem in groupitem:
                seqdf['repp'][seqdf['id'] == groupitemitem] = str(pd.Series(groupitem).values)
                seqdf['Occurrence'][seqdf['id'] == groupitemitem] = nuac[i]
            i = i +1
            
## preparing data frame for exporting as CSV file 

        seqdf.sort_values(['start'], ascending=False ,inplace=True)
        seqdf= seqdf.loc[-seqdf.iloc[:,-2].duplicated(),]
        
        seqdf.reset_index(drop=True,inplace=True)
        
        seqdf.columns = ['seqnum', 'ref1', '1', 'ref2', 'start', 'deletion', 'end', 'ref3', '4', 'ref4', 'ref4vs4', 'ref2vs4', 'ref3vs1', 'ref1vs1', 'count', 'start end', 'dwnMH', 'upMH', 'MH', 'id', 'seq to compare', 'isolates', 'Occurrence']
        
    
        seqdf_short = seqdf[['isolates','start','end','deletion','MH','Occurrence']]
        
        seqdf_short.columns = ['ID','START','END','DELETION SIZE','HOMOLOGY','OCCURRENCE']

        
        stac = []
        
        for mh in seqdf.index:
            stac.append(seqdf['start'][mh] - len(seqdf['upMH'][mh]))
        
        seqdf['start'] = stac
        
        occfig = []
        for oc in seqdf['Occurrence'].values:
            if oc == 1:
                occfig.append('')
            else:
                occfig.append("x"+ str(oc))
                
                
## Exporting data to CSV file              
                
        seqdf_short.to_csv(str(filename)+".csv")
        
## reference sequence will appear as 0 length deletion


## creating deletion spectra figures
        plt.figure(figsize=(11,len(seqsdata)/4), dpi=1000)
        filenameT = filename.split('.')[0]
        plt.text(0, len(seqdf.index)+1, str(filenameT)+ ' (n= ' + str(len(seqsdata)-1)+')' , size= 12)
        plt.axhline(y= len(seqdf.index), xmin=0, xmax=len(refseq), linewidth=10, color='red')
        
        i=0
        for loca in range(0,len(seqdf.index)-1):
            x=int(seqdf.loc[seqdf.index==loca,'start'])/int(seqdf.loc[seqdf.index[-1],'end'])
            y=int(seqdf.loc[seqdf.index==loca,'end'])/int(seqdf.loc[seqdf.index[-1],'end'])
            plt.axhline(y= len(seqdf.index), xmin=362/int(seqdf.loc[seqdf.index[-1],'end']), xmax=507/int(seqdf.loc[seqdf.index[-1],'end']), linewidth=8, color='white')
            plt.axhline(y=i, xmin=x, xmax=y, linewidth=6, color='black')

            plt.text(y+0.01, i-0.2, str(seqdf.loc[seqdf.index==loca,'MH']).split('\n')[0].split('    ')[1] + ", "+ str(int(seqdf.loc[seqdf.index==loca,'deletion'].values)) + " bp" , size= 12)
            plt.text(x-0.04, i-0.2, str(occfig[loca]) , size= 12)            
            i = i+1
            
        plt.grid(False)
        plt.axis('off')
        plt.axes().get_yaxis().set_visible(False)
        plt.savefig(str(filename)+".jpg", bbox_inches='tight', pad_inches=0)
 

