/*
A KBase module: hipmer
*/

module hipmer {
    /*
        Run assembler

        workspace_name - the name of the workspace for input/output
        read_library_name - the name of the PE read library (SE library support in the future)
        output_contig_set_name - the name of the output contigset

        extra_params - assembler specific parameters
        min_contig_length - minimum length of contigs to output, default 200

        @optional min_contig_len
        @optional extra_params
    */
    typedef structure {
        string workspace_name;
        string read_library_name;
        string output_contigset_name;

        int min_contig_len;
        list <string> extra_params;
    } AssemblyParams;

    typedef structure {
        string report_name;
        string report_ref;
    } AssemblyOutput;

    funcdef run_hipmer_hpc(AssemblyParams params) returns (AssemblyOutput output)
        authentication required;

};
