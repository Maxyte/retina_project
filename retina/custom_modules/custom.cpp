/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"
#include <bits/stdc++.h>
// #include "../../Physicell/core/PhysiCell_cell.cpp"
void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = retina_velocity_update_function;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = competence_function; // update_cell_and_death_parameters_O2_based;
	cell_defaults.functions.custom_cell_rule = NULL;
	cell_defaults.functions.contact_function = NULL;
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL;
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/*
       Cell rule definitions 
	*/

	setup_cell_rules(); 

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function;
    
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here.
    
	// initialize BioFVM
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	
	// create some of each type of cell 
	
	Cell* pC;
	
	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			position[0] = Xmin + UniformRandom()*Xrange;
			position[1] = Ymin + UniformRandom()*Yrange;
			position[2] = Zmin + UniformRandom()*Zrange;
			
			pC = create_cell( *pCD ); 
			pC->assign_position( position );
		}
	}
	std::cout << std::endl; 
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml(); 	
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

std::vector<double> L1_normalize( std::vector<double> v )
{
    double norm = 0;
    for ( int i = 0; i < v.size(); i++ )
        norm += v[i];
    if ( norm > 0.0001 )
    {
        for ( int i = 0; i < v.size(); i++ )
            v[i] = v[i] / norm;
    }
    return v;
}

void competence_function( Cell* pCell, Phenotype& phenotype, double dt ) 
{
    // ==== SETUP CONSTANTS ====
    // get current competence and normalize
    double Intermediate_C = get_single_behavior( pCell, "custom:Intermediate_competence" );
    double RGC_C = get_single_behavior( pCell, "custom:RGC_competence" );
    double MG_C = get_single_behavior( pCell, "custom:MG_competence" );
    double AC_C = get_single_behavior( pCell, "custom:AC_competence" );
    double BC_C = get_single_behavior( pCell, "custom:BC_competence" );
    double HC_C = get_single_behavior( pCell, "custom:HC_competence" );
    double Rod_C = get_single_behavior( pCell, "custom:Rod_competence" );
    double Cone_C = get_single_behavior( pCell, "custom:Cone_competence" );
    std::vector<double> competence = { Intermediate_C, RGC_C, MG_C, AC_C, BC_C, HC_C, Rod_C, Cone_C };
    competence = L1_normalize( competence );

    // TIF apices (time at which TIF is maximally expressed)
    double t = get_single_signal( pCell, "time" );
    double t0 = 0;
    double t1 = parameters.doubles["Ikzf1_temporal_marker"].value;
    double t2 = parameters.doubles["Pou2f_temporal_marker"].value;
    double t3 = parameters.doubles["Onecut_temporal_marker"].value;
    double t4 = parameters.doubles["FoxN4_temporal_marker"].value;
    double t5 = parameters.doubles["miRNA_temporal_marker"].value;
    double t6 = parameters.doubles["Nfi_temporal_marker"].value;
    double t7 = parameters.doubles["Casz1_temporal_marker"].value;
    double t8 = parameters.doubles["mature_temporal_marker"].value;
    
    // helper L1 norm function for pdf

    // ==== TEMPORAL COMPETENECE EXPECTATION ====
        // interpolate through competence vectors at TIF apices
        std::vector<double> competence_expectation;
        
        // t=0 to Ikzf1 apex
        if ( t <= t1 ) {
            std::vector<double> initial_competence_vec = {
                parameters.doubles["initial_Intermediate_competence"].value,
                parameters.doubles["initial_RGC_competence"].value,
                parameters.doubles["initial_MG_competence"].value,
                parameters.doubles["initial_AC_competence"].value,
                parameters.doubles["initial_BC_competence"].value,
                parameters.doubles["initial_HC_competence"].value,
                parameters.doubles["initial_Rod_competence"].value,
                parameters.doubles["initial_Cone_competence"].value
            };
            competence_expectation = initial_competence_vec;
        }
        // Ikzf1 to Pou2f
        if (t > t1 && t <= t2 ) {
            std::vector<double> Ikzf1_competence_vec = {
                parameters.doubles["Ikzf1_Intermediate_competence"].value,
                parameters.doubles["Ikzf1_RGC_competence"].value,
                parameters.doubles["Ikzf1_MG_competence"].value,
                parameters.doubles["Ikzf1_AC_competence"].value,
                parameters.doubles["Ikzf1_BC_competence"].value,
                parameters.doubles["Ikzf1_HC_competence"].value,
                parameters.doubles["Ikzf1_Rod_competence"].value,
                parameters.doubles["Ikzf1_Cone_competence"].value
            };
            competence_expectation = Ikzf1_competence_vec;
        }
        // Pou2f to Onecut
        if (t > t2 && t <= t3 ) {
            
            std::vector<double> Pou2f_competence_vec = {
                parameters.doubles["Pou2f_Intermediate_competence"].value,
                parameters.doubles["Pou2f_RGC_competence"].value,
                parameters.doubles["Pou2f_MG_competence"].value,
                parameters.doubles["Pou2f_AC_competence"].value,
                parameters.doubles["Pou2f_BC_competence"].value,
                parameters.doubles["Pou2f_HC_competence"].value,
                parameters.doubles["Pou2f_Rod_competence"].value,
                parameters.doubles["Pou2f_Cone_competence"].value
            };
            
            competence_expectation = Pou2f_competence_vec;
        }
        // Onecut to FoxN4
        if (t > t3 && t <= t4 ) {
            std::vector<double> Onecut_competence_vec = {
                parameters.doubles["Onecut_Intermediate_competence"].value,
                parameters.doubles["Onecut_RGC_competence"].value,
                parameters.doubles["Onecut_MG_competence"].value,
                parameters.doubles["Onecut_AC_competence"].value,
                parameters.doubles["Onecut_BC_competence"].value,
                parameters.doubles["Onecut_HC_competence"].value,
                parameters.doubles["Onecut_Rod_competence"].value,
                parameters.doubles["Onecut_Cone_competence"].value
            };
            competence_expectation = Onecut_competence_vec;
        }
        // FoxN4 to miRNA
        if (t > t4 && t <= t5 ) {
            std::vector<double> FoxN4_competence_vec = {
                parameters.doubles["FoxN4_Intermediate_competence"].value,
                parameters.doubles["FoxN4_RGC_competence"].value,
                parameters.doubles["FoxN4_MG_competence"].value,
                parameters.doubles["FoxN4_AC_competence"].value,
                parameters.doubles["FoxN4_BC_competence"].value,
                parameters.doubles["FoxN4_HC_competence"].value,
                parameters.doubles["FoxN4_Rod_competence"].value,
                parameters.doubles["FoxN4_Cone_competence"].value
            };
            competence_expectation = FoxN4_competence_vec;
        }
        // miRNA to Nfi
        if (t > t5 && t <= t6 ) {
            std::vector<double> miRNA_competence_vec = {parameters.doubles["miRNA_Intermediate_competence"].value,
                parameters.doubles["miRNA_RGC_competence"].value,
                parameters.doubles["miRNA_MG_competence"].value,
                parameters.doubles["miRNA_AC_competence"].value,
                parameters.doubles["miRNA_BC_competence"].value,
                parameters.doubles["miRNA_HC_competence"].value,
                parameters.doubles["miRNA_Rod_competence"].value,
                parameters.doubles["miRNA_Cone_competence"].value};
            competence_expectation = miRNA_competence_vec;
        }
        // Nfi to Casz1
        if (t > t6 && t <= t7 ) {
            std::vector<double> Nfi_competence_vec ={
                parameters.doubles["Nfi_Intermediate_competence"].value,
                parameters.doubles["Nfi_RGC_competence"].value,
                parameters.doubles["Nfi_MG_competence"].value,
                parameters.doubles["Nfi_AC_competence"].value,
                parameters.doubles["Nfi_BC_competence"].value,
                parameters.doubles["Nfi_HC_competence"].value,
                parameters.doubles["Nfi_Rod_competence"].value,
                parameters.doubles["Nfi_Cone_competence"].value
            };
            competence_expectation = Nfi_competence_vec;
        }
        // Casz1 to Maturity
        if ( t > t7 && t <= t8 ) {
            std::vector<double> Casz1_competence_vec = {
                parameters.doubles["Casz1_Intermediate_competence"].value,
                parameters.doubles["Casz1_RGC_competence"].value,
                parameters.doubles["Casz1_MG_competence"].value,
                parameters.doubles["Casz1_AC_competence"].value,
                parameters.doubles["Casz1_BC_competence"].value,
                parameters.doubles["Casz1_HC_competence"].value,
                parameters.doubles["Casz1_Rod_competence"].value,
                parameters.doubles["Casz1_Cone_competence"].value
            };
            
            competence_expectation = Casz1_competence_vec;
        }
        // Maturity
        if ( t > t8) 
        {
            std::vector<double> mature_competence_vec = {
                parameters.doubles["mature_Intermediate_competence"].value,
                parameters.doubles["mature_RGC_competence"].value,
                parameters.doubles["mature_MG_competence"].value,
                parameters.doubles["mature_AC_competence"].value,
                parameters.doubles["mature_BC_competence"].value,
                parameters.doubles["mature_HC_competence"].value,
                parameters.doubles["mature_Rod_competence"].value,
                parameters.doubles["mature_Cone_competence"].value
            };
            competence_expectation = mature_competence_vec;
        }
        // normalize competence expectation
        competence_expectation = L1_normalize( competence_expectation );
        
        // ==== NOTCH SIGNAL STRENGTH AND DIRECTION ====
        std::vector<double> notch_signal = {0,0,0,0,0,0,0,0};
        std::vector<Cell*> neighbors = find_nearby_interacting_cells( pCell );
        double notch_strength = 0.0;
        if ( neighbors.size() != 0 )
        {
            for ( int n = 0; n < neighbors.size(); n++ )
            {
                static std::string nb_type = neighbors[n]->type_name;
                if ( nb_type == "Intermediate" )
                    notch_signal[0] +=1;
                if ( nb_type == "RGC" )
                    notch_signal[1] +=1;
                if ( nb_type == "MG" )
                    notch_signal[2] +=1;
                if ( nb_type == "AC" )
                    notch_signal[3] +=1;
                if ( nb_type == "BC" )
                    notch_signal[4] +=1;
                if ( nb_type == "HC" )
                    notch_signal[5] +=1;
                if ( nb_type == "Rod" )
                    notch_signal[6] +=1;
                if ( nb_type == "Cone" )
                    notch_signal[7] +=1;
            }
            notch_signal = L1_normalize( notch_signal );
            auto max_nb_ratio = *max_element( notch_signal.begin(), notch_signal.end());
            notch_strength = ( max_nb_ratio - (1 / neighbors.size()) ) / ( 1 - (1 / neighbors.size()) );
        }
        // ==== COMBINE NOTCH AND TIF SIGNALS ====
    for ( int i = 0; i < competence.size(); i++ )
    {
        // equal weight at t=0 evolves linearly to TIF dominance at t_mat
        double w_TIF = (t8 + t) / (2 * t8);
        if (w_TIF > 1)
            w_TIF = 1;
        double D_TIF = w_TIF * parameters.doubles("competence_retention_factor") * (competence_expectation[i] - competence[i]);
        double w_NB = 1 - w_TIF;
        double D_NB = w_NB * notch_strength * (notch_signal[i] - competence[i]);
        // Euler's Method
        competence[i] = competence[i] + dt * ( D_TIF + D_NB );
        if (competence[i] < 0)
            competence[i] = 0;
    }
    // normalize final competence
    competence = L1_normalize( competence );
    
    set_single_behavior( pCell, "custom:Intermediate_competence", competence[0] );
    set_single_behavior( pCell, "custom:RGC_competence", competence[1] );
    set_single_behavior( pCell, "custom:MG_competence", competence[2] );
    set_single_behavior( pCell, "custom:AC_competence", competence[3] );
    set_single_behavior( pCell, "custom:BC_competence", competence[4] );
    set_single_behavior( pCell, "custom:HC_competence", competence[5] );
    set_single_behavior( pCell, "custom:Rod_competence", competence[6] );
    set_single_behavior( pCell, "custom:Cone_competence", competence[7] );
    
    // TRANSFORM TO DIFFERENTIATED NEURON BASED ON COMPETENCE
    
    // Neurogenic cell division probability as a function of time and apical proximity
    if ( pCell->custom_data["custom:flagged_to_differentiate"] == 0.0 ) 
    {
        set_single_behavior( pCell, "transform to Intermediate", 0.0 );
        set_single_behavior( pCell, "transform to RGC", 0.0 );
        set_single_behavior( pCell, "transform to MG", 0.0 );
        set_single_behavior( pCell, "transform to AC", 0.0 );
        set_single_behavior( pCell, "transform to BC", 0.0 );
        set_single_behavior( pCell, "transform to HC", 0.0 );
        set_single_behavior( pCell, "transform to Rod", 0.0 );
        set_single_behavior( pCell, "transform to Cone", 0.0 );
    }
    // prompt cell to transform
    if ( pCell->custom_data["flagged_to_migrate"] == 0 )
    {
        // halt cell
        set_single_behavior( pCell, "velocity_x", 0);
        set_single_behavior( pCell, "velocity_y", 0);
        set_single_behavior( pCell, "velocity_z", 0);
        if ( UniformRandom() < (t-t1)/(t8-t1) ) // Neurogenic division decision goes here!
        {
            if ( UniformRandom() < 0.03 )
                set_single_behavior( pCell, "custom:flagged_to_differentiate", 1 );
            if ( pCell->custom_data["flagged_to_differentiate"] == 1)
            {
                // SELECTION OF CELL TYPE FOR NEUROGENIC CELL DIVISION
                int type_index = choose_event( competence );
                std::string transformation_type_target;
                if ( type_index == 0 )
                    transformation_type_target = "Intermediate";
                if ( type_index == 1 )
                    transformation_type_target = "RGC";
                if ( type_index == 2 )
                    transformation_type_target = "MG";
                if ( type_index == 3 )
                    transformation_type_target = "AC";
                if ( type_index == 4 )
                    transformation_type_target = "BC";
                if ( type_index == 5 )
                    transformation_type_target = "HC";
                if ( type_index == 6 )
                    transformation_type_target = "Rod";
                if ( type_index == 7 )
                    transformation_type_target = "Cone";
                
                set_single_behavior( pCell, "transform to " + transformation_type_target, 1);
                set_single_behavior( pCell, "custom:flagged_to_differentiate", 0 );
                set_single_behavior( pCell, "custom:flagged_to_migrate", 1);
            }
        }
    }
    return;
}

void retina_velocity_update_function( Cell* pCell, Phenotype& phenotype, double dt )
{
    // prompt cell to migrate
    if ( pCell->custom_data["flagged_to_migrate"] == 1 )
    {
        if ( pCell->custom_data["time_migrating"] > pCell->custom_data["migration_state_length"] )
        {
            set_single_behavior( pCell, "custom:flagged_to_migrate", 0);
            set_single_behavior( pCell, "custom:time_migrating", 0);
            phenotype.mechanics.cell_cell_repulsion_strength = find_cell_definition( pCell->type_name )->phenotype.mechanics.cell_cell_repulsion_strength;
            phenotype.mechanics.cell_cell_repulsion_strength = find_cell_definition( pCell->type_name )->phenotype.mechanics.cell_cell_adhesion_strength;
        }
        else
        {
            phenotype.mechanics.cell_cell_repulsion_strength = 0;
            std::vector<double> migration_direction = normalize( pCell->position );
            double direction = 1;
            if ( pCell->type_name == "Rod" || pCell->type_name == "Cone" )
            {
                // enforce outer apical boundary
                if ( dist(pCell->position, {0,0,0}) > parameters.doubles("d_apical_mature"))
                    direction = -1;
                // migrate to PR layer
                if ( dist(pCell->position, {0,0,0}) < parameters.doubles("d_PR_layer_mature"))
                    direction = 1;
                else
                {
                    direction = UniformRandom();
                }
            }
            set_single_behavior( pCell, "velocity_x", direction * migration_direction[0]);
            set_single_behavior( pCell, "velocity_y", direction * migration_direction[1]);
            set_single_behavior( pCell, "velocity_z", direction * migration_direction[2]);
            double time_m = get_single_behavior( pCell, "custom:time_migrating") + dt;
            set_single_behavior( pCell, "time_migrating", time_m );
        }
    }
}

void dividing_phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{
    return;
}

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; }

//        // REFERENCE FOR RPC COMPETENCE AT EACH TIF APEX
//        std::vector<double> initial_competence_vec = {
//            parameters.doubles["initial_Intermediate_competence"].value,
//            parameters.doubles["initial_RGC_competence"].value,
//            parameters.doubles["initial_MG_competence"].value,
//            parameters.doubles["initial_AC_competence"].value,
//            parameters.doubles["initial_BC_competence"].value,
//            parameters.doubles["initial_HC_competence"].value,
//            parameters.doubles["initial_Rod_competence"].value,
//            parameters.doubles["initial_Cone_competence"].value
//        };
//        std::vector<double> Ikzf1_competence_vec = {
//            parameters.doubles["Ikzf1_Intermediate_competence"].value,
//            parameters.doubles["Ikzf1_RGC_competence"].value,
//            parameters.doubles["Ikzf1_MG_competence"].value,
//            parameters.doubles["Ikzf1_AC_competence"].value,
//            parameters.doubles["Ikzf1_BC_competence"].value,
//            parameters.doubles["Ikzf1_HC_competence"].value,
//            parameters.doubles["Ikzf1_Rod_competence"].value,
//            parameters.doubles["Ikzf1_Cone_competence"].value
//        };
//
//        std::vector<double> Pou2f_competence_vec = {
//            parameters.doubles["Pou2f_Intermediate_competence"].value,
//            parameters.doubles["Pou2f_RGC_competence"].value,
//            parameters.doubles["Pou2f_MG_competence"].value,
//            parameters.doubles["Pou2f_AC_competence"].value,
//            parameters.doubles["Pou2f_BC_competence"].value,
//            parameters.doubles["Pou2f_HC_competence"].value,
//            parameters.doubles["Pou2f_Rod_competence"].value,
//            parameters.doubles["Pou2f_Cone_competence"].value
//        };
//
//        std::vector<double> Onecut_competence_vec = {
//            parameters.doubles["Onecut_Intermediate_competence"].value,
//            parameters.doubles["Onecut_RGC_competence"].value,
//            parameters.doubles["Onecut_MG_competence"].value,
//            parameters.doubles["Onecut_AC_competence"].value,
//            parameters.doubles["Onecut_BC_competence"].value,
//            parameters.doubles["Onecut_HC_competence"].value,
//            parameters.doubles["Onecut_Rod_competence"].value,
//            parameters.doubles["Onecut_Cone_competence"].value
//        };
//
//        std::vector<double> FoxN4_competence_vec = {
//            parameters.doubles["FoxN4_Intermediate_competence"].value,
//            parameters.doubles["FoxN4_RGC_competence"].value,
//            parameters.doubles["FoxN4_MG_competence"].value,
//            parameters.doubles["FoxN4_AC_competence"].value,
//            parameters.doubles["FoxN4_BC_competence"].value,
//            parameters.doubles["FoxN4_HC_competence"].value,
//            parameters.doubles["FoxN4_Rod_competence"].value,
//            parameters.doubles["FoxN4_Cone_competence"].value
//        };
//
//        std::vector<double> miRNA_competence_vec = {
//            parameters.doubles["miRNA_Intermediate_competence"].value,
//            parameters.doubles["miRNA_RGC_competence"].value,
//            parameters.doubles["miRNA_MG_competence"].value,
//            parameters.doubles["miRNA_AC_competence"].value,
//            parameters.doubles["miRNA_BC_competence"].value,
//            parameters.doubles["miRNA_HC_competence"].value,
//            parameters.doubles["miRNA_Rod_competence"].value,
//            parameters.doubles["miRNA_Cone_competence"].value
//        };
//
//        std::vector<double> Nfi_competence_vec ={
//            parameters.doubles["Nfi_Intermediate_competence"].value,
//            parameters.doubles["Nfi_RGC_competence"].value,
//            parameters.doubles["Nfi_MG_competence"].value,
//            parameters.doubles["Nfi_AC_competence"].value,
//            parameters.doubles["Nfi_BC_competence"].value,
//            parameters.doubles["Nfi_HC_competence"].value,
//            parameters.doubles["Nfi_Rod_competence"].value,
//            parameters.doubles["Nfi_Cone_competence"].value
//        };
//
//        std::vector<double> Casz1_competence_vec = {
//            parameters.doubles["Casz1_Intermediate_competence"].value,
//            parameters.doubles["Casz1_RGC_competence"].value,
//            parameters.doubles["Casz1_MG_competence"].value,
//            parameters.doubles["Casz1_AC_competence"].value,
//            parameters.doubles["Casz1_BC_competence"].value,
//            parameters.doubles["Casz1_HC_competence"].value,
//            parameters.doubles["Casz1_Rod_competence"].value,
//            parameters.doubles["Casz1_Cone_competence"].value
//        };
//        std::vector<double> mature_competence_vec = {
//            parameters.doubles["mature_Intermediate_competence"].value,
//            parameters.doubles["mature_RGC_competence"].value,
//            parameters.doubles["mature_MG_competence"].value,
//            parameters.doubles["mature_AC_competence"].value,
//            parameters.doubles["mature_BC_competence"].value,
//            parameters.doubles["mature_HC_competence"].value,
//            parameters.doubles["mature_Rod_competence"].value,
//            parameters.doubles["mature_Cone_competence"].value
//        };

