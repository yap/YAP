#ifndef __examples_bat_find_pdl_file
#define __examples_bat_find_pdl_file

/// \return path to default .pdl file
inline const std::string find_pdl_file()
{ return (::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") + "/data" : "./data") + "/evt.pdl"; }

#endif
