#include "neb.h"
#include "neb_fortran_interface.h"

using namespace pele;

namespace {
	NEB *g_neb = nullptr;
}

void neb_setup(potential_callback_t *potential, void *userdata)
{
	if (g_neb)
		throw std::logic_error("neb already initialized");
	g_neb = new NEB(nullptr);
}

void neb_cleanup()
{
	if (!g_neb)
		throw std::logic_error("trying to cleanup neb which was not initialized");
	delete g_neb;
	g_neb = nullptr;
}

void neb_start()
{
	g_neb->start();
}

bool neb_step()
{
	return g_neb->step();
}



