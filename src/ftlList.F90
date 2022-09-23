! Instantiate the module from the ftlList template
#define FTL_TEMPLATE_TYPE integer
#define FTL_TEMPLATE_TYPE_NAME Int
#define FTL_INSTANTIATE_TEMPLATE
#include "ftlList.F90_template"

! Make fpm recognize the instantiated module(needs to be fixed in fpm)
#if 0
module ftlListIntModule
end module ftlListIntModule
#endif