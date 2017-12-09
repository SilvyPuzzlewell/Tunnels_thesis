
#ifndef PATH_DUPLUCATION_CHECK
#define PATH_DUPLICATION_CHECK

#include "world.h"
#include <memory>

int is_tunnel_duplicated(shared_ptr<Path> checked_path, std::vector<shared_ptr<Path>> existing_paths);

#endif 