##################################################################################################
#
# TracSynth.py - GW quality trace element analysis & Monte Carlo modeling setup
#
# (1) clean up (account for non-detects, etc.)
# (2) limit to useful analytes (e.g., sufficient detections)
# (3) compute and list correlation matrix
# (4) compute covariance matrix and posit correlated random parameter values (e.g., for PHREEQC)
#
##################################################################################################

from numpy import *
from pandas import *
import matplotlib.pyplot as plt
import subprocess
import os


def Tally(results, sorbed, analytes):
    # tally adsorbed masses per components in analytes list
    for colName in list(results):
        results.rename(columns={colName: colName.replace('m_', '')}, inplace=True)
    sorbedResults = results[sorbed].copy()
    for i, chem in enumerate(analytes):
        results[chem + '-sorbed'] = 0.
        for sorb in sorbed:
            if chem in sorb: results[chem + '-sorbed'] = results[chem + '-sorbed'] + results[sorb]
    # drop individual adsorbed species from results
    results.drop(sorbed, axis=1, inplace=True)    
    return results


def ReadSorbed():
    # read text file with list of adsorbed species
    lineInput = []
    sorbed = []        
    inputFile = open('surface_species.txt','r')
    for line in inputFile: lineInput.append(line.split()[0])
    inputFile.close()
    for i in range(len(lineInput)): sorbed.append(lineInput[i])
    print('Read surface species listing.')     
    return sorbed


def ReadWQ(numSampleMin):
    # read water quality data and return data frame
    wq = read_csv('water_quality.csv')
    analytes = list(wq)
    analytes.remove('Well')
    analytes.remove('Date')    
    wq['Date'] = to_datetime(wq['Date'])
    ndCount = zeros(len(analytes), int)    
    for j, chem in enumerate(analytes):
        conc = list(wq[chem].astype('str'))
        for i in range(len(conc)):
            # delete R, J, and U flags (comments in analytical data)
            if conc[i].find('R')!=-1:
                conc[i] = conc[i].replace('R', '')
                conc[i] = conc[i].strip()
            if conc[i].find('J-')!=-1:
                conc[i] = conc[i].replace('J-', '')
                conc[i] = conc[i].strip()                
            if conc[i].find('J')!=-1:
                conc[i] = conc[i].replace('J', '')
                conc[i] = conc[i].strip()
            if conc[i].find('U')!=-1:
                conc[i] = conc[i].replace('U', '')
                conc[i] = conc[i].strip()                
            # convert U to 0.5 x DL, or remove sample, depending on modeDL
            if conc[i].find('<')!=-1:
                conc[i] = conc[i].replace('<', '')
                conc[i] = conc[i].strip()
                conc[i] = 0.5 * float(conc[i])
                ndCount[j] += 1
        wq[chem] = array(conc).astype('float')
        if chem != 'pH': wq[chem] = log10(wq[chem])
    wq.dropna(axis=0, how='all', subset=analytes, inplace=True)
    # remove analytes with too many non-detects from consideration    
    samples = wq.count()
    samples.drop(labels=['Well', 'Date'], inplace=True)
    samples = samples-ndCount
    samples = samples[samples>numSampleMin]
    validAnalytes = list(samples.keys())
    names = validAnalytes
    header = ['Well', 'Date']
    header.extend(validAnalytes)
    print('Read and processed water quality data.')
    return wq[header], names
    

def CorrelPlots(corrData, corrSynth, names, analytes, indx):
    # generate correlation vs correlation plots (synthetic versus observed)
    for i, chem in enumerate(analytes):
        plt.figure(indx)
        x = array(corrData[chem])
        y = array(corrSynth[chem])
        plt.scatter(x, y, s=15, facecolors='blue', edgecolors='blue')
        for j, name in enumerate(names):
            plt.annotate(name, (x[j]+0.02, y[j]))
        plt.title(chem)
        plt.xlabel('Correlation (observations)')
        plt.ylabel('Correlation (synthetic)')
        plt.show()


def Scatter(data, synth, xSpecies, ySpeciesSet):
    # generate selected scatter plots
    for i, chem in enumerate(ySpeciesSet):
        plt.figure(i)
        plt.scatter(data[xSpecies], data[chem], s=25, facecolors='black', edgecolors='black', label = 'Data')    
        plt.scatter(synth[xSpecies], synth[chem], s=3, facecolors='none', edgecolors='red', label = 'Synthetic')
        plt.xscale('log')
        plt.xlabel(xSpecies + ' (mol/L)')
        if chem != 'pH':
            plt.ylabel(chem + ' (mol/L)')
            plt.yscale('log')
        else:
            plt.ylabel(chem)
        plt.legend(loc=2)
        plt.show()


def WriteInput(synthetic, names, equilParam, equilMin, analytes, special, sorbed, phases):

    # write out PHREEQC input file
    N = len(synthetic)
    phrqFile = open('phrqInput.txt','w')
    phrqFile.writelines(['TITLE GW Speciation Model', '\n'])

    # solution spread keyword block
    phrqFile.writelines(['', '\n'])    
    phrqFile.writelines(['SOLUTION_SPREAD', '\n'])    
    phrqFile.writelines(['\t', '-units', '\t', 'ug/l', '\n'])
    header = '\t'
    equils = '\t'
    for i, chem in enumerate(names):
        header += chem + '\t'        
        for j, param in enumerate(equilParam):
            if chem == param: equils += equilMin[j] + '  0' + '\t'  # fix concentrations via phase equilibrium
            else: equils += '\t'
    phrqFile.writelines([header, '\n']) 
    phrqFile.writelines([equils, '\n']) 
    for i in range(N):
        row = '\t'
        for j in range(len(names)): row += str(synthetic[names[j]].iloc[i]) + '\t'
        phrqFile.writelines([row, '\n'])
        
    # surface blocks    
    for i in range(N): 
        phrqFile.writelines(['', '\n'])
        phrqFile.writelines(['SURFACE ' + str(i+1), '\n'])
        phrqFile.writelines(['\t', '-equilibrate with solution', '\t', str(i+1), '\n'])
        phrqFile.writelines(['\t', 'Hfo_w', '\t', '0.005', '\t', '600', '\t', '1', '\n'])
        phrqFile.writelines(['\t', 'Hfo_s', '\t', '0.00005', '\n'])    

    # equilibrate with phases across all solutions; use run_cells to distribute
    phrqFile.writelines(['', '\n']) 
    phrqFile.writelines(['EQUILIBRIUM_PHASES', '\t', str(1) + '-' + str(N), '\n'])
    for phase in phases:
        phrqFile.writelines([phase, '\t', str(0.), '\t', str(0.), '\n'])
    phrqFile.writelines(['', '\n'])
    phrqFile.writelines(['RUN_CELLS', '\n'])
    phrqFile.writelines(['\t', '-cells', '\t', str(1) + '-' + str(N), '\n'])

    # selected output block
    phrqFile.writelines(['', '\n'])    
    phrqFile.writelines(['SELECTED_OUTPUT', '\n'])
    phrqFile.writelines(['\t', '-file', '\t', 'selected_output.txt', '\n'])
    phrqFile.writelines(['\t', '-reset', '\t', 'false', '\n'])
    phrqFile.writelines(['\t', '-state', '\t', 'true', '\n'])
    phrqFile.writelines(['\t', '-solution', '\t', 'true', '\n'])
    phrqFile.writelines(['\t', '-pH', '\t', 'true', '\n'])
    phrqFile.writelines(['\t', '-pe', '\t', 'true', '\n'])
    phrqFile.writelines(['\t', '-percent_error', '\t', 'true', '\n'])
    analytesString = ['\t', '-totals']
    for chem in special:
        analytesString.append('\t')
        analytesString.append(chem)        
    for chem in analytes:
        analytesString.append('\t')
        analytesString.append(chem)
    analytesString.append('\n')
    phrqFile.writelines(analytesString)
    sorbedString =  ['\t', '-molalities']
    for sorb in sorbed:
        sorbedString.append('\t')
        sorbedString.append(sorb)
    sorbedString.append('\n')
    phrqFile.writelines(sorbedString)
	# note that tracked saturation indices do not impact model and may not be needed; hard-wired here for expediency
    phrqFile.writelines(['\t', '-saturation_indices', '\t', 'Calcite  CO2(g)  Hydroxylapatite  FCO3Apatite', '\n'])
    phrqFile.writelines(['', '\n'])
    phrqFile.writelines(['END', '\n'])    
    phrqFile.close()

    
### main script ###    

def TraceSynth(N):

    numSampleMin = 5                    # minimum number of detections to include analyte
    wq, names = ReadWQ(numSampleMin)    # process water quality observational data  
    
    # species that will be processed as output
    components = read_csv('components.csv')     
    analytes = list(components['Analyte'])
    MW = list(components['MW'])
    compsDict = dict(zip(analytes, MW))
    
    # additonal parameters to be included in output
    sorbed = ReadSorbed()
    sorbedElements = ['P', 'S', 'As', 'Ca', 'Cr', 'Pb', 'Mg', 'V'] 		# track total amounts adsorbed
    special = ['As(3)', 'As(5)']
    diffTrack = ['As', 'P', 'pH'] 			# track updated concentrations following mineral equilibration reactions

    # phases for equilibrium phase blocks
    phases = ['Hydroxylapatite']

    # note correlations between parameters in data
    corrWQ = DataFrame(wq.corr())
    corrWQ.to_csv('data_corr.csv') 

    # exclude analytes that are not sufficiently co-reported with other analaytes
    exclude = ['Mn', 'TPH-g', 'TPH-r', 'TPH-d'] 	# hard-wired list of analytes to be dropped from analysis
    wq.drop(exclude, axis=1, inplace=True)    
    for x in exclude: names.remove(x)       # names = set of useable parameter observations
    
    # generate synthetic data (from log ug/L concentrations)
    mean0 = wq.mean().values
    cov0 = wq.cov().values
    r = random.multivariate_normal(mean0, cov0, N)
    synthetic =DataFrame(r)
    synthetic.columns = names
    synthetic['pe'] = 4.0
    names.append('pe')

    # convert synthetic data to linear units (ug/L)
    for chem in names:
        if (chem != 'pH') and (chem != 'pe'):
            wq[chem] = 10**wq[chem]
            synthetic[chem] = 10**synthetic[chem]

    # remove posited solutions that exceed concentration caps
    clipComp = ['Cl']                   # components list and associated maximum allowed concentration
    clipMax = [3000.*1000.]
    for i, chem in enumerate(clipComp):
        synthetic = synthetic[synthetic[chem]<=clipMax[i]].copy()

    # run modeled solutions, surfaces, and equilibrium phases in a single batch
    WriteInput(synthetic, names, ['pe'], ['Ferrihydrite'], analytes, special, sorbed, phases)                                  # write PHREEQC input file
    subprocess.call('PHREEQC phrqInput.txt output.txt minteq.v4.dat', shell=True)           # run PHREEQC
    
    # process model output; drop results with very poor charge imbalances (e.g., greater than 40%) 
    selected = read_csv('selected_output.txt', delim_whitespace=True) 
    selected.drop_duplicates(inplace=True)    
    results0 = selected[(selected['state']=='i_surf') & (abs(selected['pct_err'])<=40)].copy()
    results0.drop('state', axis=1, inplace=True)    
    results0 = Tally(results0, sorbed, analytes)
    resultsM = selected[(selected['state']=='react') & (abs(selected['pct_err'])<=40)].copy()
    resultsM.drop('state', axis=1, inplace=True)      
    results0.to_csv('synthetic_results.csv', index=False)
    
    # process output for graphical representations
    analytes.append('pH')
    syntheticFiltered = results0[analytes]   # filtered to include specific analytes of interest (with molecular weights)
    for chem in analytes:       # convert water quality data to mol/L to facilitate comparison
        if chem!='pH': wq[chem] = 1e-6 * wq[chem] / compsDict.get(chem)

    # show example observed & synthetic data scatter plots
    yShifts = zeros(len(analytes), float)
    Scatter(wq, syntheticFiltered, 'As', analytes)

    # write screened observational data to file
    wq.to_csv('wq_set.csv', index=False)

    # convert back to log space
    synthSubset = results0[analytes].copy()
    for chem in analytes:
        if chem != 'pH':
            wq[chem] = log10(wq[chem])
            synthSubset[chem] = log10(synthSubset[chem])

    # display correlation-to-correlation plots
    corrData = DataFrame(wq.corr())
    corrSynth = DataFrame(synthSubset.corr())
    indx = len(analytes)    
    CorrelPlots(corrData, corrSynth, analytes, ['Fe'], indx)
    
    # merge in post-prepipitation results; display selected parameter distributions to graphs
    keep = ['soln']
    keep.extend(diffTrack)
    resultsM = resultsM[keep]
    for chem in diffTrack: resultsM.rename(columns={chem: chem + '_m'}, inplace=True)
    resultsAll = results0.merge(resultsM, on='soln', how='inner')
    resultsAll.to_csv('results_all.csv', index=False)
    resultsAll['diffP'] = resultsAll['P'] - resultsAll['P_m']
    resultsAll['rAs'] = (resultsAll['As'] - resultsAll['As_m']) / resultsAll['As']  
    resultsAll['rP'] = (resultsAll['P'] - resultsAll['P_m']) / resultsAll['P']    

	# scatter plot comparing initial arsenic concentrations versus those corresponding to P removal
    plt.figure(indx)
    resultsAll.plot.scatter(x='As', y='As_m', s=resultsAll['diffP'] * 100000)
    plt.xlabel('Arsenic (observed)')
    plt.ylabel('Arsenic (corrected)')
    plt.xscale('log')     
    plt.yscale('log')    
    plt.show()

	# histogram of arsenic concentration reductions (from adsorption)
    plt.figure(indx+2)
    hist = resultsAll['rAs'].hist(bins=20)
    plt.title('Arsenic Concentration Changes')
    plt.xlabel('Concentration Reduction (fraction)')
    plt.ylabel('N')    
    plt.show()

	# histogram of phosphate concentration reductions (from hydroxylapatite precipitation)	
    plt.figure(indx+3)
    hist = resultsAll['rP'].hist(bins=20)
    plt.title('Phosphate Concentration Changes')    
    plt.xlabel('Concentration Reduction (fraction)')
    plt.ylabel('N')     
    plt.show()

    print('N = ', len(resultsAll))
    print('Done.')
    

### run script ###

N = 750 		# requested size of synthetic population; may drop because of implied charge imbalances
TraceSynth(N)        
    
    