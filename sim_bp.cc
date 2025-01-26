#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstdint>
#include "sim_bp.h"
#include <math.h>
class Branch_predictor{
    public:
    int no_of_counters;
    int* bp_array;
    int no_of_bits_gb;
    int predictions = 0;
    int mispredictions = 0 ;
    uint32_t global_branch_history_n;

    Branch_predictor(uint32_t m, int n  ) : no_of_counters(m), no_of_bits_gb(n) { 
        int power = pow(2,no_of_counters);
        bp_array = new int[power];
        global_branch_history_n = 0;
        intialize(power);

    }

    virtual void intialize(int power){
        for(int i = 0; i < power; ++i){
            bp_array[i] = 2; // initialized to 2 (“weakly taken”)
        }
    }

    virtual ~Branch_predictor(){
        delete[] bp_array;
    }

    virtual int predictor(uint32_t b_index){
        //printf("prediction done");

        return (bp_array[b_index] >= 2) ? 1 : 0;
        //printf("prediction done");        
    } 

    virtual int get_index(uint32_t b_index) {


        return b_index;

    }

    virtual bool actual_outcome_changes(int act_out,int value_pre , uint32_t b_index){

        if((act_out == 1) && (bp_array[b_index] < 3)){
            bp_array[b_index]++;
            if(bp_array[b_index] > 3){
                bp_array[b_index] = 3;
            }
            //printf("bp++");
        }
        else if((act_out == 0) && (bp_array[b_index] > 0)){
            bp_array[b_index]--;
            if(bp_array[b_index] < 0){
                bp_array[b_index] = 0;
            }
            //printf("bp--");
        }
        else if((bp_array[b_index] > 3) || (bp_array[b_index] < 0))
        {
            printf("2_bit_counter out of bound");
        }

        predictions++;
        if(act_out != value_pre){
            mispredictions++;
            return false;            
        }
        else {
            return true;
        }
    }
    virtual bool found_in_array(int act_out, int value_pre){
        //predictions++;
        if(act_out != value_pre){
            //mispredictions++;
            return false;            
        }
        else {
            return true;
        }
    }

    uint32_t m_bits_addr(uint32_t addr){

        return (addr >> 2) & ((1 << no_of_counters) - 1);

    }

    virtual bool read_addr(uint32_t addr, int act_out,Branch_predictor*bimodel,Branch_predictor*gshare){
        uint32_t index = m_bits_addr(addr);
        index = get_index(index);
        int value_pre = predictor(index);
        bool found = actual_outcome_changes(act_out, value_pre ,index);
        global_branch_history(act_out);

        return found;

        //printf("done");

    }
    virtual void global_branch_history(int act_out){
        if(no_of_bits_gb > 0){

        global_branch_history_n = (global_branch_history_n >> 1) | (act_out << (no_of_bits_gb-1));
        //printf(".");
        
        }
    }
    void display(){
            for (int j = 0; j < pow(2,no_of_counters); ++j) {
    
                printf(" %d\t%d\n",j ,bp_array[j]);
    
        }
    
    }
    void printMeasurements(){
    double ratio = (((double)mispredictions/(double)predictions)*100.0);
    printf("OUTPUT\n");
    printf(" number of predictions:    %d\n",predictions);
    printf(" number of mispredictions: %d\n",mispredictions);
    printf(" misprediction rate:       %.2f%%\n",ratio);
}
};

class bimodel : public Branch_predictor{
public:
    bimodel(uint32_t m, int n) : Branch_predictor(m,n) {}
        
};

class gshare : public Branch_predictor{
public:
    gshare(uint32_t m, int n) : Branch_predictor(m,n) {}


    int get_index(uint32_t b_index) override{

    uint32_t m_n = no_of_counters - no_of_bits_gb;
    uint32_t mask = ((1U << no_of_bits_gb) - 1) << m_n;
    
    uint32_t shifted_gb = global_branch_history_n << m_n;
    uint32_t result = b_index ^ (shifted_gb & mask);

    return result;

    }

};

class hybrid : public Branch_predictor{
public:
    hybrid(uint32_t m, int n) : Branch_predictor(m,n) {
        int power = pow(2,no_of_counters);
        bp_array = new int[power];
        global_branch_history_n = 0;
        intialize(power);
    }

    void intialize(int power) override{
        for(int i = 0; i < power; ++i){
            bp_array[i] = 1; // initialized to 2 (“weakly taken”)
        }
    }
    bool read_addr(uint32_t addr, int act_out,Branch_predictor*bimodel,Branch_predictor*gshare) override{
        uint32_t index = m_bits_addr(addr);
        //index = get_index(index);
        //int value_pre = predictor(index);
        predictions++;
        bool found_in_bi = bimodel->found_in_array(act_out,bimodel->predictor(bimodel->get_index(bimodel->m_bits_addr(addr))));
        bool found_in_gs = gshare->found_in_array(act_out,gshare->predictor(gshare->get_index(gshare->m_bits_addr(addr))));
        bool found = false;


        if(bp_array[index] >= 2){
            found = gshare->read_addr(addr,act_out,nullptr,nullptr);
        }
        else if(bp_array[index] <= 1){
            found = bimodel->read_addr(addr,act_out,nullptr,nullptr);
            gshare->global_branch_history(act_out);
        }
        if(!found){
            mispredictions++;
        }

            if(!found_in_bi && found_in_gs){
                bp_array[index]++;
                if(bp_array[index] > 3){
                bp_array[index]=3;
            }
            }
            else if(found_in_bi && !found_in_gs){
                bp_array[index]--;
                if(bp_array[index] < 0){
                bp_array[index]=0;
            }
            }

        return false;
        //printf("done");

    }


        
};



/*  argc holds the number of command line arguments
    argv[] holds the commands themselves

    Example:-
    sim bimodal 6 gcc_trace.txt
    argc = 4
    argv[0] = "sim"
    argv[1] = "bimodal"
    argv[2] = "6"
    ... and so on
*/
int main (int argc, char* argv[])
{
    FILE *FP;               // File handler
    char *trace_file;       // Variable that holds trace file name;
    bp_params params;       // look at sim_bp.h header file for the the definition of struct bp_params
    char outcome;           // Variable holds branch outcome
    unsigned long int addr; // Variable holds the address read from input file

    
    if (!(argc == 4 || argc == 5 || argc == 7))
    {
        printf("Error: Wrong number of inputs:%d\n", argc-1);
        exit(EXIT_FAILURE);
    }
    
    params.bp_name  = argv[1];

   
    
    // strtoul() converts char* to unsigned long. It is included in <stdlib.h>
    if(strcmp(params.bp_name, "bimodal") == 0)              // Bimodal
    {
        if(argc != 4)
        {
            printf("Error: %s wrong number of inputs:%d\n", params.bp_name, argc-1);
            exit(EXIT_FAILURE);
        }
        params.M2       = strtoul(argv[2], NULL, 10);
        trace_file      = argv[3];
        printf("COMMAND\n%s %s %lu %s\n", argv[0], params.bp_name, params.M2, trace_file);
        

    }
    else if(strcmp(params.bp_name, "gshare") == 0)          // Gshare
    {
        if(argc != 5)
        {
            printf("Error: %s wrong number of inputs:%d\n", params.bp_name, argc-1);
            exit(EXIT_FAILURE);
        }
        params.M1       = strtoul(argv[2], NULL, 10);
        params.N        = strtoul(argv[3], NULL, 10);
        trace_file      = argv[4];
        printf("COMMAND\n%s %s %lu %lu %s\n", argv[0], params.bp_name, params.M1, params.N, trace_file);

    }
    else if(strcmp(params.bp_name, "hybrid") == 0)          // Hybrid
    {
        if(argc != 7)
        {
            printf("Error: %s wrong number of inputs:%d\n", params.bp_name, argc-1);
            exit(EXIT_FAILURE);
        }
        params.K        = strtoul(argv[2], NULL, 10);
        params.M1       = strtoul(argv[3], NULL, 10);
        params.N        = strtoul(argv[4], NULL, 10);
        params.M2       = strtoul(argv[5], NULL, 10);
        trace_file      = argv[6];
        printf("COMMAND\n%s %s %lu %lu %lu %lu %s\n", argv[0], params.bp_name, params.K, params.M1, params.N, params.M2, trace_file);

    }
    else
    {
        printf("Error: Wrong branch predictor name:%s\n", params.bp_name);
        exit(EXIT_FAILURE);
    }
    
    // Open trace_file in read mode
    FP = fopen(trace_file, "r");
    if(FP == NULL)
    {
        // Throw error and exit if fopen() failed
        printf("Error: Unable to open file %s\n", trace_file);
        exit(EXIT_FAILURE);
    }
    bimodel bi(params.M2, 0);
    gshare  gs(params.M1, params.N);
    hybrid  hy(params.K,0);
    bool found;
    

    char str[2];
    while(fscanf(FP, "%lx %s", &addr, str) != EOF)
    {
    if(strcmp(params.bp_name, "bimodal") == 0){
        outcome = str[0];
        if (outcome == 't'){
            //printf("%lx %s\n", addr, "t");           // Print and test if file is read correctly
            found = bi.read_addr(addr,1,nullptr,nullptr);
        }
        else if (outcome == 'n'){
            //printf("%lx %s\n", addr, "n");          // Print and test if file is read correctly
            found = bi.read_addr(addr,0,nullptr,nullptr);
        }
    }
    else if (strcmp(params.bp_name, "gshare") == 0){
        outcome = str[0];
        if (outcome == 't'){
            //printf("%lx %s\n", addr, "t");           // Print and test if file is read correctly
            found = gs.read_addr(addr,1,nullptr,nullptr);
        }
        else if (outcome == 'n'){
            //printf("%lx %s\n", addr, "n");          // Print and test if file is read correctly
            found = gs.read_addr(addr,0,nullptr,nullptr);
        }
    }
    else if (strcmp(params.bp_name, "hybrid") == 0){
        outcome = str[0];
        if (outcome == 't'){
            //printf("%lx %s\n", addr, "t");           // Print and test if file is read correctly
            found = hy.read_addr(addr,1,&bi,&gs);
        }
        else if (outcome == 'n'){
            //printf("%lx %s\n", addr, "n");          // Print and test if file is read correctly
            found = hy.read_addr(addr,0,&bi,&gs);
        }
    }
    }

    if(strcmp(params.bp_name, "bimodal") == 0){
        bi.printMeasurements();
        printf("FINAL BIMODAL CONTENTS\n");
        bi.display();
        
    }
    else if (strcmp(params.bp_name, "gshare") == 0){
        gs.printMeasurements();
        printf("FINAL GSHARE CONTENTS\n");
        gs.display();

    }
    else if (strcmp(params.bp_name, "hybrid") == 0){
        hy.printMeasurements();
        printf("FINAL CHOOSER CONTENTS\n");
        hy.display();
        printf("FINAL GSHARE CONTENTS\n");
        gs.display();
        printf("FINAL BIMODAL CONTENTS\n");
        bi.display();
    }

    return 0;
}
