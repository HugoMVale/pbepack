! Instantiate 'ftlListIntModule' from the ftlList template
#define FTL_TEMPLATE_TYPE integer
#define FTL_TEMPLATE_TYPE_NAME Int
#define FTL_INSTANTIATE_TEMPLATE
#include "ftlList.F90_template"

! Instantiate 'ftlListDoubleModule' from the ftlList template
#define FTL_TEMPLATE_TYPE real(kind=8)
#define FTL_TEMPLATE_TYPE_NAME Double
#define FTL_INSTANTIATE_TEMPLATE
#include "ftlList.F90_template"

! Make fpm recognize the instantiated modules (needs to be fixed in fpm)
#if 0
module ftlListIntModule
end module ftlListIntModule
module ftlListDoubleModule
end module ftlListDoubleModule
#endif