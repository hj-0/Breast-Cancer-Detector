#define _CRTDBG_MAP_ALLOC
#include "DecisionTree.h"
#include "FileReaderWriter.h"
#include <cstdlib>
#include <crtdbg.h>

using namespace std;

int main(int argc, char* argv[])
{
	// check for any memory leaks
#if defined(_DEBUG)
	int dbgFlags = _CrtSetDbgFlag(_CRTDBG_REPORT_FLAG);
	dbgFlags |= _CRTDBG_CHECK_ALWAYS_DF;
	dbgFlags |= _CRTDBG_DELAY_FREE_MEM_DF;
	dbgFlags |= _CRTDBG_LEAK_CHECK_DF;
	_CrtSetDbgFlag(dbgFlags);

	cout << "\t\t\t\t\tBreast Cancer Diagnosis with Binary Decision Tree" << endl;
	cout << "\t\t\t\t\t*************************************************" << endl;
#endif

	vector<DecisionTree::WisconsinNode*> v_attributes;
	auto* tree = new DecisionTree();
	FileReaderWriter f;

	f.read_csv(tree, v_attributes);

	int no_of_malignant = 0;
	int no_of_benign = 0;
	for (auto& v : v_attributes) {
		if (v->class_ == 2) {
			++no_of_benign;
		}
		else if (v->class_ == 4) {
			++no_of_malignant;
		}
	}

	cout << "\n\nTotal Patients processed: " << no_of_benign + no_of_malignant << endl;
	cout << "Total Benign: " << no_of_benign << endl;
	cout << "Total Malignant: " << no_of_malignant << endl;

	f.write_csv(v_attributes);

	tree->clear_vector(v_attributes); // empty the vector along with the nodes
	delete tree; // delete a ptr

	return 0;
}