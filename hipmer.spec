/*
A KBase module: hipmer
*/

module hipmer {
    /*
        Run assembler

        workspace_name - the name of the workspace for input/output
        read_library_name - the name of the PE read library (SE library support in the future)
        output_contig_set_name - the name of the output contigset

        mer_size -
        is_diploid - is gnome diploid

    */
    typedef structure {
        string workspace_name;
        string read_library_name;
        string output_contigset_name;

        int mer_size;
        int is_diploid;
        float dynamic_min_depth;
        int min_depth_cutoff;
        float gap_close_rpt_depth_ratio;
        int assm_scaff_len_cutoff;
    } AssemblyParams;

    typedef structure {
        string report_name;
        string report_ref;
    } AssemblyOutput;

    funcdef run_hipmer_hpc(UnspecifiedObject params) returns (AssemblyOutput output)
        authentication required;

};
