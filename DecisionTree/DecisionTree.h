#pragma once
#include <vector>

class DecisionTree {
public:
	struct WisconsinNode {
		// all attributes initialized to 0
		int id_ = 0, clump_thickness_ = 0, uniformity_of_cell_size_ = 0, uniformity_of_cell_shape_ = 0, marginal_adhesion_ = 0,
			single_epithelial_cell_size_ = 0, bare_nuclei_ = 0, bland_chromatin_ = 0, normal_nucleoli_ = 0, mitoses_ = 0, class_ = 0;

		// roots
		WisconsinNode* wisconsin_root_ = nullptr;   // wisconsin root node initialized to nullptr
		WisconsinNode* wisconsin_left_ = nullptr;   // wisconsin left child node initialized to nullptr
		WisconsinNode* wisconsin_right_ = nullptr;  // wisconsin right child node initialized to nullptr

		//Default ctor
		WisconsinNode() = default;

		// Copy ctor
		WisconsinNode(WisconsinNode * attribute) {
			id_ = attribute->id_;
			clump_thickness_ = attribute->clump_thickness_;
			uniformity_of_cell_size_ = attribute->uniformity_of_cell_size_;
			uniformity_of_cell_shape_ = attribute->uniformity_of_cell_shape_;
			marginal_adhesion_ = attribute->marginal_adhesion_;
			single_epithelial_cell_size_ = attribute->single_epithelial_cell_size_;
			bare_nuclei_ = attribute->bare_nuclei_;
			bland_chromatin_ = attribute->bland_chromatin_;
			normal_nucleoli_ = attribute->normal_nucleoli_;
			mitoses_ = attribute->mitoses_;
			class_ = attribute->class_;
		}
	};

private:
	WisconsinNode* root_; // this is the first node in our binary decision tree
public:
	//default ctor with root_ variable initialized to nullptr
	DecisionTree() : root_(nullptr) {}

	// operations
	WisconsinNode* construct_decision_tree(WisconsinNode* node);
	static void delete_nodes(WisconsinNode* node, std::vector<WisconsinNode*>& v_attributes);
	static void clear_vector(std::vector<WisconsinNode*>& v_attributes);

};
//delete all nodes in the binary decision tree recursively
inline void DecisionTree::delete_nodes(WisconsinNode* node, std::vector<WisconsinNode*>& v_attributes) {
	if (node == nullptr) {
		return;
	}
	if (node->wisconsin_left_ != nullptr) {
		delete_nodes(node->wisconsin_left_, v_attributes);
	}
	if (node->wisconsin_right_ != nullptr) {
		delete_nodes(node->wisconsin_right_, v_attributes);
	}

	delete node;
}

//clear the vector along with the nodes stored recursively
inline void DecisionTree::clear_vector(std::vector<WisconsinNode*>& v_attributes) {
	for (size_t i = 0; i < v_attributes.size(); ++i) {
		delete_nodes(v_attributes[i], v_attributes);
	}
}

// Binary decision tree which determines which types of breast cancer tumors a patient is diagnosed with
inline DecisionTree::WisconsinNode* DecisionTree::construct_decision_tree(WisconsinNode* node) {
	this->root_ = node;
	node->wisconsin_root_ = node;

	if (node->uniformity_of_cell_size_ <= 2) {
		auto* new_node_a = new WisconsinNode(node);
		node->wisconsin_left_ = new_node_a;
		node = node->wisconsin_left_;

		if (node->bare_nuclei_ <= 3) {
			auto* new_node_b = new WisconsinNode(node);
			node->wisconsin_left_ = new_node_b;
			node = node->wisconsin_left_;
			node->class_ = 2;
			this->root_->class_ = 2;  // assign class_ value to root_ of 'DecisionTree' obj
			return this->root_; // return class_ value to 'DecisionTree' obj for later storing it into the vector in main scope
		}
		if (node->bare_nuclei_ > 3) {
			auto* new_node_c = new WisconsinNode(node);
			node->wisconsin_right_ = new_node_c;
			node = node->wisconsin_right_;

			if (node->clump_thickness_ <= 3) {
				auto* new_node_d = new WisconsinNode(node);
				node->wisconsin_left_ = new_node_d;
				node = node->wisconsin_left_;
				node->class_ = 2;
				this->root_->class_ = 2;
				return this->root_;
			}
			if (node->clump_thickness_ > 3) {
				auto* new_node_e = new WisconsinNode(node);
				node->wisconsin_right_ = new_node_e;
				node = node->wisconsin_right_;

				if (node->bland_chromatin_ <= 2) {
					auto* new_node_f = new WisconsinNode(node);
					node->wisconsin_left_ = new_node_f;
					node = node->wisconsin_left_;

					if (node->marginal_adhesion_ <= 3) {
						auto* new_node_g = new WisconsinNode(node);
						node->wisconsin_left_ = new_node_g;
						node = node->wisconsin_left_;
						node->class_ = 4;
						this->root_->class_ = 4;
						return this->root_;
					}
					if (node->marginal_adhesion_ > 3) {
						auto* new_node_h = new WisconsinNode(node);
						node->wisconsin_right_ = new_node_h;
						node = node->wisconsin_right_;
						node->class_ = 2;
						this->root_->class_ = 2;
						return this->root_;
					}
				}
				else if (node->bland_chromatin_ > 2) {
					auto* new_node_i = new WisconsinNode(node);
					node->wisconsin_right_ = new_node_i;
					node = node->wisconsin_right_;
					node->class_ = 4;
					this->root_->class_ = 4;
					return this->root_;
				}
			}
		}
	}
	else if (node->uniformity_of_cell_size_ > 2) {
		auto* new_node_a = new WisconsinNode(node);
		node->wisconsin_right_ = new_node_a;
		node = node->wisconsin_right_;

		if (node->uniformity_of_cell_shape_ <= 2) {
			auto* new_node_b = new WisconsinNode(node);
			node->wisconsin_left_ = new_node_b;
			node = node->wisconsin_left_;

			if (node->clump_thickness_ <= 5) {
				auto* new_node_c = new WisconsinNode(node);
				node->wisconsin_left_ = new_node_c;
				node = node->wisconsin_left_;
				node->class_ = 2;
				this->root_->class_ = 2;
				return this->root_;
			}
			if (node->clump_thickness_ > 5) {
				auto* new_node_d = new WisconsinNode(node);
				node->wisconsin_right_ = new_node_d;
				node = node->wisconsin_right_;
				node->class_ = 4;
				this->root_->class_ = 4;
				return this->root_;
			}
		}
		else if (node->uniformity_of_cell_shape_ > 2) {
			auto* new_node_e = new WisconsinNode(node);
			node->wisconsin_right_ = new_node_e;
			node = node->wisconsin_right_;

			if (node->uniformity_of_cell_size_ <= 4) {
				auto* new_node_f = new WisconsinNode(node);
				node->wisconsin_left_ = new_node_f;
				node = node->wisconsin_left_;

				if (node->bare_nuclei_ <= 2) {
					auto* new_node_g = new WisconsinNode(node);
					node->wisconsin_left_ = new_node_g;
					node = node->wisconsin_left_;

					if (node->marginal_adhesion_ <= 3) {
						auto* new_node_h = new WisconsinNode(node);
						node->wisconsin_left_ = new_node_h;
						node = node->wisconsin_left_;
						node->class_ = 2;
						this->root_->class_ = 2;
						return this->root_;
					}
					if (node->marginal_adhesion_ > 3) {
						auto* new_node_i = new WisconsinNode(node);
						node->wisconsin_right_ = new_node_i;
						node = node->wisconsin_right_;
						node->class_ = 4;
						this->root_->class_ = 4;
						return this->root_;
					}
				}
				else if (node->bare_nuclei_ > 2) {
					auto* new_node_j = new WisconsinNode(node);
					node->wisconsin_right_ = new_node_j;
					node = node->wisconsin_right_;

					if (node->clump_thickness_ <= 6) {
						auto* new_node_k = new WisconsinNode(node);
						node->wisconsin_left_ = new_node_k;
						node = node->wisconsin_left_;

						if (node->uniformity_of_cell_size_ <= 3) {
							auto* new_node_l = new WisconsinNode(node);
							node->wisconsin_left_ = new_node_l;
							node = node->wisconsin_left_;
							node->class_ = 4;
							this->root_->class_ = 4;
							return this->root_;
						}
						if (node->uniformity_of_cell_size_ > 3) {
							auto* new_node_m = new WisconsinNode(node);
							node->wisconsin_right_ = new_node_m;
							node = node->wisconsin_right_;

							if (node->marginal_adhesion_ <= 5) {
								auto* new_node_n = new WisconsinNode(node);
								node->wisconsin_left_ = new_node_n;
								node = node->wisconsin_left_;
								node->class_ = 2;
								this->root_->class_ = 2;
								return this->root_;
							}
							if (node->marginal_adhesion_ > 5) {
								auto* new_node_o = new WisconsinNode(node);
								node->wisconsin_right_ = new_node_o;
								node = node->wisconsin_right_;
								node->class_ = 4;
								this->root_->class_ = 4;
								return this->root_;
							}
						}
					}
					else if (node->clump_thickness_ > 6) {
						auto* new_node_p = new WisconsinNode(node);
						node->wisconsin_right_ = new_node_p;
						node = node->wisconsin_right_;
						node->class_ = 4;
						this->root_->class_ = 4;
						return this->root_;
					}
				}
			}
			else if (node->uniformity_of_cell_size_ > 4) {
				auto* new_node_q = new WisconsinNode(node);
				node->wisconsin_right_ = new_node_q;
				node = node->wisconsin_right_;
				node->class_ = 4;
				this->root_->class_ = 4;
				return this->root_;
			}
		}
	}
	return this->root_;
}







