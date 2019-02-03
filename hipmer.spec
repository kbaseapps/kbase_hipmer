/*
A KBase module: hipmer
*/

module hipmer {

    typedef structure {
        string report_name;
        string report_ref;
    } AssemblyOutput;

    funcdef run_hipmer_hpc(UnspecifiedObject params) returns (AssemblyOutput output)
        authentication required;

};
