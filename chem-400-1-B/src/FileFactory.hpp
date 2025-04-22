#ifndef FILE_FACTORY
#define FILE_FACTORY

#include <sstream>
#include <fstream>

#include "Frame.hpp"
#include "ConfigHandler.hpp"
#include "Mesh.hpp"
#include "Chemotherapy.hpp"
#include "../src_backup/PharmacologicModel.hpp"

class FileFactory {
public:
	static void makeFrameFile(Frame* frame, ConfigHandler* config);
	static void makeFile(Frame* frame, ConfigHandler* config, Mesh* mesh, Chemotherapy* chemotherapy, OutputMode mode, int chemDrug = -1);
	static void makePKFile(PharmacologicModel* pharmacologic, ConfigHandler* config);
	static void makePDFile(PharmacologicModel* pharmacologic, ConfigHandler* config);
};

#endif /* end of include guard: FILE_FACTORY */
