
mutable struct PerturbSensitivityModel
    model::LinearSensitivityModel
    sp::AbstractArray #Real power voltage magnitude sensitivity matrix
    sq::AbstractArray #Reactive power " "
    spq::AbstractArray #Interleaved matrix
    ΔV::AbstractArray #Estimated Voltage Deviations
    function PerturbSensitivityModel(system_model::System)
        model = LinearSensitivityModel(system_model)
        sp = calc_dv_dp_matrix(system_model)
        sq = cal_dv_dq_matrix(system_model)
        spq = calc_sens_mat_perturb(system_model)
        return new(model,sp,sq,spq,ΔV) 
    end
end

function calc_sens_mat_perturb(sensors,injections,system_model::System)
    """Makes a perturb-and-observe sensitivity matrix for a model"""
    M = size(sensors)
    L = size(injections)
    S = zeros(M,L)
    for (sensor,injection) in zip(sensors,injections)
        1+2
        #voltages_inj = injection_voltages()
    end
end

# function calc_dv_dp_matrix(sensors,injections,system_model::System)
#     """Form the sensitivity matrix for a distribution model.

#     Keyword arguments:
#     opendssDir -- The folder the .dss file is located in
#     opendssFile -- The .dss file name
#     circuit_name -- the "short name" for the circuit e.g. "ieee13"
#     """
#     M = size(sensors)
#     L = size(injections)
#     S = zeros(M,L)

#     os.chdir(opendssDir)
#     dss.run_command('Compile "'+opendssFile+'"')
#     dss.run_command('Set Controlmode = STATIC')
#     dss.run_command('solve')
#     ######
#     # Options for truncating the number of candidate nodes: 
#     #   node_range - lets you manually set the range of nodes
#     #   phase_mode - takes a brute_force approach and cuts off the feeder
#     #   node_select - lets you manually select candidate nodes.
    
        
#     if phase_mode == '3ph':
#         buses = dss.Circuit.AllBusNames()[3:]
#         #Get the base case voltages
#         voltages_basecase = dss.Circuit.AllBusMagPu()[3:]
    
    
#     elseif phase_mode == '1ph' :
#         buses = dss.Circuit.AllNodeNames()[9:]
#         voltages_basecase = dss.Circuit.AllBusMagPu()[9:]
    
#     else:
#         buses = truncate_vector(dss.Circuit.AllNodeNames(),circuit_name)
#         voltages_basecase = truncate_vector(dss.Circuit.AllBusMagPu(),circuit_name) 
#         ####
#         ## Remove the electricall identical bus and also the low side of the transformer. 
#         ## A new method is needed to calculate the columns that are not necessary so that we can 
#         ## Generalize to a bigger feeder case.
#         ## 
#         ## ADD MORE LATER
#         ##
#         ## ROUTINE 1: IEEE13. SINGLE PHASE Electrically identical buses were established previously using MATLAB.
#         ##if('ieee13' in circuit_name):   
#         ##    buses = [x for x in buses[9:] if "634" not in x if "692" not in x]
            
        
#     end 
#     #Preprocess the injection labels
    
#     injections = []
#     j=0
#     k=0
#     for ix in range(2*len(buses)):
#         if(np.remainder(ix,2)!=0):
#             injections.append([buses[k],buses[k]+' Q'])
#             k+=1
#         else:
#             injections.append([buses[j],buses[j]+' P'])
#             j+=1
    
#     #if('ieee13' in circuit_name):
#     #    injections = [x for x in injections if "634" not in x if "692" not in x]
    
#     #Iteratively solve the injections and retrieve the voltages 
#     S_PQ = np.zeros((len(buses),len(injections)))
    
#     for col_idx, injection in enumerate(injections):
        
#         if('P' in injection[1]):
#             voltages_inj = injection_voltages(opendssDir,'injection',injection[0],'1','-100','0',opendssFile=opendssFile)
#             voltages_inj = truncate_vector(voltages_inj,circuit_name,idx_select=idx_select)
            
#         if('Q' in injection[1]):
#             voltages_inj= injection_voltages(opendssDir,'injection',injection[0],'1','0','-100',opendssFile=opendssFile)
#             voltages_inj = truncate_vector(voltages_inj,circuit_name,idx_select=idx_select)
            
#         S_PQ[:,col_idx] = (np.asarray(voltages_inj)-np.asarray(voltages_basecase))/100
        
#     injection_labels = [injection[1] for injection in injections]
    
#     return S_PQ,pd.DataFrame(data=S_PQ,index=buses,columns=injection_labels),injection_labels,voltages_basecase

# end

