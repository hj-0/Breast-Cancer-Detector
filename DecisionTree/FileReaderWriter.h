#pragma once
#include "DecisionTree.h"
#include <iostream>
#include <string>
#include <fstream>

using namespace std;

class FileReaderWriter {
public:
	// operations
	void read_csv(DecisionTree* tree, vector<DecisionTree::WisconsinNode*>& v_attributes) const; // read from csv file
	static void write_csv(vector<DecisionTree::WisconsinNode*>& v_attributes); // write to csv file
	
};

inline void FileReaderWriter::read_csv(DecisionTree* tree, vector<DecisionTree::WisconsinNode*>& v_attributes) const {

	ifstream csv("unformatted_data_v1.0.0.csv", ios::in);
	if (!csv) {
		cout << "\nError: failed to open 'unformatted_data_v1.0.0.csv'" << endl;
	}

	string id, clump_thickness, uniformity_of_cell_size, uniformity_of_cell_shape, marginal_adhesion, single_epithelial_cell_size,
		bare_nuclei, bland_chromatin, normal_nucleoli, mitoses, class_column;

	cout << "reading 'unformatted_data_v1.0.0.csv'...";

	/*int invalid_patients = 0;*/
	while (csv.good()) {
		while (getline(csv, id, ',')) {
			getline(csv, clump_thickness, ',');
			getline(csv, uniformity_of_cell_size, ',');
			getline(csv, uniformity_of_cell_shape, ',');
			getline(csv, marginal_adhesion, ',');
			getline(csv, single_epithelial_cell_size, ',');
			getline(csv, bare_nuclei, ',');
			getline(csv, bland_chromatin, ',');
			getline(csv, normal_nucleoli, ',');
			getline(csv, mitoses, ',');
			getline(csv, class_column, '\n');

			DecisionTree::WisconsinNode* attributes = new DecisionTree::WisconsinNode();

			// this will skip patients which has missing data
			// but for some unknown reason I am getting memory leaks from this so I decided to replace "?" with "1" to get rid of the memory leaks in the next step
			/*if (bare_nuclei == "?") {
				++invalid_patients;
				break;
			}*/

			if(bare_nuclei == "?") { // if exception is thrown, replace "?" with "1"
				bare_nuclei = "1";
			}

			attributes->id_ = stoi(id);
			attributes->clump_thickness_ = stoi(clump_thickness);
			attributes->uniformity_of_cell_size_ = stoi(uniformity_of_cell_size);
			attributes->uniformity_of_cell_shape_ = stoi(uniformity_of_cell_shape);
			attributes->marginal_adhesion_ = stoi(marginal_adhesion);
			attributes->single_epithelial_cell_size_ = stoi(single_epithelial_cell_size);
			attributes->bare_nuclei_ = stoi(bare_nuclei);
			attributes->bland_chromatin_ = stoi(bland_chromatin);
			attributes->normal_nucleoli_ = stoi(normal_nucleoli);
			attributes->mitoses_ = stoi(mitoses);
			attributes->class_ = stoi(class_column);

			attributes = tree->construct_decision_tree(attributes);

			v_attributes.push_back(attributes);
		}
	}

	/*cout << "\nTotal Patients with missing data: " << invalid_patients << endl;*/
	csv.close();
}

inline void FileReaderWriter::write_csv(vector<DecisionTree::WisconsinNode*>& v_attributes)  {

	ofstream csv("results.csv", ios::out);

	if(!csv) {
		cout << "\nError: failed to write to 'results.csv'" << endl;
	}

	cout << "\nwriting to 'results.csv'..." << endl;
	cout << "please wait..." << endl;
	for(auto& v : v_attributes) {	
		csv << v->clump_thickness_ << ", " << v->uniformity_of_cell_size_ << "," << v->uniformity_of_cell_shape_ << ","
			<< v->marginal_adhesion_ << "," << v->single_epithelial_cell_size_ << "," << v->bare_nuclei_ << "," << v->bland_chromatin_ << ","
			<< v->normal_nucleoli_ << "," << v->mitoses_ << "," << v->class_ << endl;
	}
	
	csv.close();
	cout << "done." << endl;
}


