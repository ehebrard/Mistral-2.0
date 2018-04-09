#ifndef MTREE_H
#define    MTREE_H

#include <cmath>
#include <vector>

#include "XCSP3utils.h"

using namespace std;
namespace XCSP3Mistral {



    //-------------------------------------

    typedef enum {
        NT_CONSTANT, NT_VARIABLE, NT_NEG, NT_ABS, NT_ADD, NT_SUB, NT_MULT, NT_DIV, NT_MOD, NT_MIN, NT_MAX, NT_LE, NT_DIST, NT_LT, NT_GE, NT_GT, NT_EQ,
        NT_NE, NT_NOT, NT_AND, NT_OR, NT_XOR, NT_IF, NT_IFF, NT_IMP, NT_IN, NT_SET
    } NodeType;

    class Node {

    public:
        NodeType type;

        Node(NodeType n) : type(n) {}

        virtual ~Node() {}
    };


    //-------------------------------------

    class NodeConstant : public Node {
    public:
        int val;


        NodeConstant(int v) : Node(NT_CONSTANT), val(v) {}

    };

    //-------------------------------------

    class NodeVariable : public Node {
    public:
        string var;


        NodeVariable(string v) : Node(NT_VARIABLE), var(v) {}
    };
		
    // //-------------------------------------
    //
    // class NodeSet : public Node {
    // public:
    //     vector<int> values;
    //
    //
    //     NodeSet(vector<int>& vals) : Node(NT_SET), values(vals) {}
    // };


    //-------------------------------------

    class NodeOperator : public Node {
    public:
        vector<Node *> args;
        int maxParameters;
        NodeOperator(NodeType n, int maxParams = 2) : Node(n), maxParameters(maxParams) {}


        void addParameter(Node *p) {
            args.push_back(p);
            if(maxParameters != -1 && (int)(args.size()) > maxParameters)
                throw runtime_error("Too many parameters in expression");
        }

        virtual ~NodeOperator() {
            for(Node *n : args)
                delete n;
        }
    };



    static int min3(int v1, int v2, int v3) {
        if(v1 == -1) v1 = INT_MAX;
        if(v2 == -1) v2 = INT_MAX;
        if(v3 == -1) v3 = INT_MAX;

        if(v1 < v2) {
            if(v1 < v3) return v1;
            return v3;
        }
        if(v2 < v3) return v2;
        return v3;
    }


    class Tree {
    protected:
        std::string expr;

    public:
        Node *root;


        Tree(std::string e) : expr(e) {
            root = fromStringToTree(expr);
        }

        void dispose() {
            delete root;
        }


        Node *fromStringToTree(std::string current) {
            current = XCSP3Core::trim(current);
            //    int pos = current.find('(');
            vector<NodeOperator *> stack;
            vector<Node *> params;
            while(true) {
                int posOpenParenthesis = current.find('(');
                int posCloseParenthesis = current.find(')');
                int posComma = current.find(',');
                int nb = min3(posCloseParenthesis, posComma, posOpenParenthesis);


                string currentElement = current.substr(0, nb);

                if(currentElement != "" && nb != posOpenParenthesis)
                    createBasicParameter(currentElement, stack, params);


                if(nb == posCloseParenthesis)
                    closeOperator(stack, params);


                if(nb == posOpenParenthesis)
                    createOperator(currentElement, stack, params);

                current = current.substr(nb + 1);
                if(current == "") break;
            }
            assert(params.size() == 1);
            assert(stack.size() == 0);
            return params.back();
        }

        void createBasicParameter(string currentElement,vector<NodeOperator*> &stack,vector<Node*> &params) {
            try {
                int nb = stoi(currentElement);
                //cout << "Number " << nb << endl;
                params.push_back(new NodeConstant(nb));
            } catch (invalid_argument e) {
							
                params.push_back(new NodeVariable(currentElement));
                /*
                 * bool alreadyInsideList = false;
                for (int i = 0; i < listOfVariables.size(); i++)
                    if (listOfVariables[i] == v) {
                        alreadyInsideList = true;
                        break;
                    }
                if (alreadyInsideList == false) listOfVariables.push(v);
                 */
							
            }
        }

        void createOperator(string currentElement, vector<NodeOperator *> &stack, vector<Node *> &params) {
            NodeOperator *tmp = NULL;
            if(currentElement == "neg") tmp = new NodeOperator(NT_NEG, 1);
            if(currentElement == "abs") tmp = new NodeOperator(NT_ABS, 1);

            if(currentElement == "add") tmp = new NodeOperator(NT_ADD, -1);
            if(currentElement == "sub") tmp = new NodeOperator(NT_SUB);
            if(currentElement == "mul") tmp = new NodeOperator(NT_MULT, -1);
            if(currentElement == "div") tmp = new NodeOperator(NT_DIV);
            if(currentElement == "mod") tmp = new NodeOperator(NT_MOD);

            //if (currentElement == "sqr") tmp = new NodeSquare(SQR);
            //if (currentElement == "pow") tmp = new NodePow(POW);

            if(currentElement == "min") tmp = new NodeOperator(NT_MIN, -1);
            if(currentElement == "max") tmp = new NodeOperator(NT_MAX, -1);
            if(currentElement == "dist") tmp = new NodeOperator(NT_DIST);
            if(currentElement == "le") tmp = new NodeOperator(NT_LE);
            if(currentElement == "lt") tmp = new NodeOperator(NT_LT);
            if(currentElement == "ge") tmp = new NodeOperator(NT_GE);
            if(currentElement == "gt") tmp = new NodeOperator(NT_GT);

            if(currentElement == "ne") tmp = new NodeOperator(NT_NE);
            if(currentElement == "eq") tmp = new NodeOperator(NT_EQ, -1);

            if(currentElement == "not") tmp = new NodeOperator(NT_NOT, 1);
            if(currentElement == "and") tmp = new NodeOperator(NT_AND, -1);
            if(currentElement == "or") tmp = new NodeOperator(NT_OR, -1);
            if(currentElement == "xor") tmp = new NodeOperator(NT_XOR);
            if(currentElement == "imp") tmp = new NodeOperator(NT_IMP);
            if(currentElement == "if") tmp = new NodeOperator(NT_IF,3);
            if(currentElement == "iff") tmp = new NodeOperator(NT_IFF);
						
						if(currentElement == "in") tmp = new NodeOperator(NT_IN);
						if(currentElement == "set") tmp = new NodeOperator(NT_SET, -1);

            //if (currentElement == "in")  throw runtime_error("IN not allowed in expression");
            //if (currentElement == "set") throw runtime_error("SET not allowed in expression");

            if(tmp == NULL)
                throw runtime_error("Intension constraint. Unknown operator: " + currentElement);
            stack.push_back(tmp);
            params.push_back(NULL); // delemitor
        }


        void closeOperator(vector<NodeOperator *> &stack, vector<Node *> &params) {
            NodeOperator *tmp = stack.back();

            int startParams = params.size() - 1;
            while(params[startParams] != NULL)
                startParams--;
            startParams++;
            int nbP = 0;
            for(size_t i = startParams ; i < params.size() ; i++, nbP++)
                tmp->addParameter(params[i]);
            stack.pop_back();
            params.resize(params.size() - nbP);
            assert(params.back() == NULL);
            params.pop_back();
            params.push_back(tmp);
            //cout << "Close operator " << "nbP = " << nbP << endl;
        }


    };
}
#endif	/* MTREE_H */
