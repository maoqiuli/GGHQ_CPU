#define PLATG   true
#define RECORD_PATH_DISTANCE  false
#define PATH_LEN 100

#define FMTINT  true
#define FMTUINT  false

#if FMTUINT
typedef uint8_t DTSET;
#else
typedef int8_t  DTSET;
#endif

#if FMTINT
typedef int     DTRES;
#else
typedef float   DTSET;
typedef float   DTRES;
#endif