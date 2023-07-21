// For data inputs // these numbers are arbitrary for the most part, but ensure enough memory
#define MAX_SUMMARY_COLS 121
#define MAX_SUMMARY_ROWS 2001
#define PARAMFILE_MAXROWS 2001
#define PARAMFILE_MAXCOLS 101
#define DATAFILE_MAXROWS 2000001
#define DATAFILE_MAXCOLS 101
#define CONFIGFILE_MAXROWS 10
#define CONFIGFILE_MAXCOLS 30

// For staging the model
#define STAGE_ID_NONE 0 // use this mode to skip all of the stage code and just run based on param sheet settings (a "normal" run, the default config)
#define STAGE_ID_HIST_OPT 1
#define STAGE_ID_HIST_STRESS 2
#define STAGE_ID_FUT_OPT 3
#define STAGE_ID_FUT_STRESS 4
#define STAGE_ID_FUT_STRESS_NOACCLIM 5