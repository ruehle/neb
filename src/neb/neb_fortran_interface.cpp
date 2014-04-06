#include "neb.h"
#include "neb_fortran_interface.h"

using namespace pele;

namespace {
	NEB *g_neb = nullptr;

}

void neb_setup(potential_callback_t *potential, void *userdata)
{
	if (g_neb)
		neb_cleanup();
	g_neb = new NEB(nullptr);
}

void neb_cleanup()
{
	if (g_neb)
		delete g_neb;
	g_neb = nullptr;
}



