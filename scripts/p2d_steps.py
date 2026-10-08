"""Canonical module registry shared by orchestration and status reporting."""
STEPS=[
('M01','01_init_project.py','M01_init.log'),
('M02','02_normalize_inputs.py','M02_normalization.log'),
('M03','03_map_pockets.py','M03_pocket_mapping.log'),
('M04','04_build_job_matrix.py','M04_job_builder.log'),
('M05A','05_engine_preflight.py','M05_engine_preflight.log'),
('M05B','05_prepare_inputs.py','M05_prepare_inputs.log'),
('M05C','05c_run_screening.py','M05C_screening.log'),
('M06','06_parse_screening.py','M06_parsing.log'),
('M07','07_rank_candidates.py','M07_ranking.log'),
('M08','08_chemoinformatics.py','M08_chemoinformatics.log'),
('M09','09_candidate_gate.py','M09_candidate_gate.log'),
('M10A','10a_build_dynamics_jobs.py','M10A_dynamics_jobs.log'),
('M10B','10b_dynamics_preflight.py','M10B_dynamics_preflight.log'),
('M10C','10c_extract_selected_pose.py','M10C_pose_extraction.log'),
('M11','11_parameterize_system.py','M11_parameterization.log'),
('M12','12_solvate_minimize.py','M12_minimization.log'),
('M13','13_equilibrate.py','M13_equilibration.log'),
('M14','14_production_md.py','M14_production_md.log'),
('M15','15_preprocess_trajectory.py','M15_trajectory_preprocessing.log'),
('M16','16_analyze_trajectory.py','M16_trajectory_analysis.log'),
('M17','17_binding_energetics.py','M17_energetics.log'),
('M18','18_essential_dynamics.py','M18_pca.log'),
('M19','19_free_energy_landscape.py','M19_landscape.log'),
('M20','20_interaction_persistence.py','M20_interactions.log'),
('M21','21_reproducibility_bundle.py','M21_bundle.log')]

def log_name(module,name,profile):
    return name.replace('.log','_fixture.log') if profile=='fixture' and module in [s[0] for s in STEPS[13:]] else name
