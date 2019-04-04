#/bin/python
import os, sys, argparse
import numpy as np

def collect_scan_results(Jparam, scanfolder):
    '''(str, str) -> dict

    Searches scanfolder for the results of Jparam scan and returns the results as a dictionary
    
    '''
    results = {}
    results[Jparam.lower()] = {}
   
    allfolder= os.listdir(scanfolder)
    for mainfolder in allfolder:

        if 'jedi-scan-' in mainfolder and os.path.isdir(scanfolder+'/'+mainfolder):
            mor = mainfolder.split('-')
            parameter = mor[-2]
            value = mor[-1]
            if parameter.lower() != Jparam.lower():
                print('ERROR: Folder %s does not contain results for a values scan of parameter %s' % (scanfolder, Jparam))
                print('ERROR: The assessment chekcpoint expected %s but found %s' % (Jparam, parameter))
                print('STOP')
                sys.exit()
            else:
                parameter = Jparam.lower()
                #initialize results dictionay for parameter 
                

            results_folder = os.listdir(scanfolder+'/'+mainfolder)
            
            for queryfolder in results_folder:
                if os.path.isdir(scanfolder+'/'+mainfolder+'/'+queryfolder):
                    # Identify which probe was used in the JEDI calculation from the name of the diferent folders
                    mor = queryfolder.split('-')
                    probe = mor[-1]
                    try:
                        for queryfile in os.listdir(scanfolder+'/'+mainfolder+'/'+queryfolder):
                            if 'jedi_stats' in queryfile:
                                stats_file = scanfolder+'/'+mainfolder+'/'+queryfolder+'/'+queryfile
                                jedistats = read_stats_file(stats_file)
                                # Identify the state of the binding site (open/collapsed or other)
                                if stats_file.endswith('initial.dat'):
                                    BSstate='open'
                                elif stats_file.endswith('14800.dat'):
                                    BSstate='collapsed'
                                else:
                                    print('WARNING: Unknown state of the protein Binding Site')
                                    BSstate='Unknown'
                                # Add the current parameter combination if it does not exists in the results dictionary

                                # First identify open, closed or Unknown
                                if BSstate not in results[parameter]:
                                    results[parameter][BSstate] = {}
                                # Second key is the used probe:
                                if probe not in results[parameter][BSstate]:
                                    results[parameter][BSstate][probe] = {}
                                # Third key is the value of the parameter
                                if value not in results[parameter][BSstate][probe]:
                                    results[parameter][BSstate][probe][value] = {}
                                #Now populate the dictionary
                                results[parameter][BSstate][probe][value]['JEDI'] = jedistats[0]
                                results[parameter][BSstate][probe][value]['Va'] = jedistats[1]
                                results[parameter][BSstate][probe][value]['Ha'] = jedistats[2]
                    except FileNotFoundError:
                        print("ERROR: File %s not found")
                        sys.exit()
    return results

def read_stats_file(jedistats):
    '''(str) -> tupple

    Reads JEDI results from the first line of jedistats
    
    '''
    jedi = 0.0
    Va = 0.0
    Ha = 0.0
    try:
        with open(jedistats) as f:
            while True:
                line = f.readline()
                if line == '':
                    break
                elif '#' in line:
                    continue
                else:
                    mor = line.split()
                    jedi = mor[1]
                    Va = mor[3]
                    Ha = mor[4]
                    return (jedi,Va,Ha)
    except:
        raise


def reformat_results(Jparam, scan_results, threshold, output='JEDI', diff=['open', 'collapsed']):
    '''(str, dict, float, ['JEDI' | 'Va' | 'Ha'], [reference_state, query_state]) -> float
    
    '''
    ref_state = diff[0]
    query_state = diff[1]

    all_plots = {}
    parameter = Jparam.lower()

    if ref_state not in scan_results[parameter]:
        print("ERROR: Results for the reference state %s not found within the gathered results" % ref_state )
        print("ERROR: The available options are")
        for kk in scan_results[parameter]:
            print('\t%s' % kk)
        print('STOP')
        sys.exit()

    if query_state not in scan_results[parameter]:
        print("ERROR: Results for the query state %s not found within the gathered results" % query_state)
        print("ERROR: The available options are")
        for kk in scan_results[parameter]:
            print('\t%s' % kk)
        print('STOP')
        sys.exit()

    ref_results = scan_results[parameter][ref_state] 
    query_results = scan_results[parameter][query_state]

    for probe in ref_results:
        if probe not in query_results:
            print('INFO: Values for probe %s were not found in results gathered for query state %s' % (probe, query_state))
            continue   
        result_key = "%s-%s" % (parameter, probe)
        paramvaluesR = []
        outvaluesR = []
        paramvaluesQ= []
        outvaluesQ= []

        for k, v in sorted(ref_results[probe].items()):
            paramvaluesR.append(float(k))
            outvaluesR.append(float(v[output]))
            if k in query_results[probe]:
                paramvaluesQ.append(float(k))
                outvaluesQ.append(float(query_results[probe][k][output]))
            else:
                print("INFO: Values for the combination of parameter: %s, value:%s, probe: %s were not found for query state: %s" % (Jparam, k, probe, query_state))
        
        all_plots[result_key] = {}
        all_plots[result_key][ref_state] = (np.array(paramvaluesR), np.array(outvaluesR))
        all_plots[result_key][query_state] = (np.array(paramvaluesQ), np.array(outvaluesQ))

    return all_plots

def check_differences (formated_results, sysid,threshold, diff=['open', 'collapsed']):
    '''(dict, str)

    '''
    try:
        ref_sys = formated_results[sysid][diff[0]]
        quer_sys = formated_results[sysid][diff[1]]
    except KeyError:
        print("ERROR: It was not possible to calculate differences between state %s and state %s for system %s" % (diff[0], diff[1], sysid) )
        raise

    try:
        temp = ref_sys[1] - quer_sys[1]
        if max(temp) > threshold:
            print("JEDI RESULT: System %s displayed a maximum difference > %.2f between %s and %s" % (sysid, threshold, diff[0], diff[1]))
            print(temp)
        else:
            print("System %s DID NOT displayed differences > %.2f between %s and %s " % (sysid, threshold, diff[0], diff[1]))
            print(temp)
    except:
        raise
    

    


def plot_diff_vs_paramvalue(scan_results):
    '''(dict) -> NoneType

    '''
    




if __name__ == '__main__':
    #Arguments to incorporate in a parser
    parser = argparse.ArgumentParser(description="Analyses the results of a scan of JEDI parameters")
    parser.add_argument('-foldername', metavar='-f', type=str, required= True, help="Folder containing results of a scan of JEDI parameters")
    parser.add_argument('-parameter', metavar='-p', type=str, required=True, help="Name of the parameter to gather results for (case insensitive)")
    parser.add_argument('-threshold', metavar='-th', type=float, required=True, help="The minimum significant difference expected in the JEDI result")
    
    args=parser.parse_args()
    foldername = args.foldername
    parameter = args.parameter
    threshold = args.threshold
    #End of parser

    scan_results = collect_scan_results(parameter,foldername)
    result_summary = reformat_results(parameter, scan_results, threshold,)
    for system in result_summary:
        check_differences(result_summary,system,threshold)
    
    #TBC: plot results for systems displaying significant differences

    
