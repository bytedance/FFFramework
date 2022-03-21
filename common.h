#if defined(FFFRAMEWORK_EXPORTS) 
#   define FFFRAMEWORK_API   __declspec(dllexport)
#else 
#   define FFFRAMEWORK_API   __declspec(dllimport)
#endif