//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: inherit.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include "scatmech.h"
#include "inherit.h"
#include "askuser.h"
#include <sstream>

#include "filmtran.h"
#include "scattabl.h"
#include "dielfunc.h"

using namespace std;


namespace SCATMECH {


    bool Inheritance::show_model_mode=false;
    int Model::quiet=0;

    //
    // stringlist takes a string and parses it into a list of
    // strings without any whitespace...
    //
    void stringlist(const string& s,deque<string>& result)
    {
        result.clear();
        result.push_back("");

        for (string::const_iterator ps=s.begin(); ps!=s.end(); ++ps) {
            if (isspace(*ps)) {
                if (!result.back().empty()) result.push_back("");
            } else {
                result.back() += *ps;
            }
        }
        if (result.back().empty()) result.pop_back();
    }

    //
    // Wrap_String takes a string and formats it into lines with an indentation.
    //
    void Wrap_String(const string& in,string& out,int indent,int width)
    {
        deque<string> ls;
        stringlist(in,ls);

        string indenter="\n";
        for (int i=0; i<indent; ++i) indenter += space;
        int wide = width-indent;

        int currentwidth = 0;
        for (deque<string>::const_iterator pls=ls.begin(); pls!=ls.end(); ++pls) {
            if (((int)pls->size()+currentwidth+1)>=wide) {
                out += indenter;
                currentwidth = pls->size()+1;
            } else {
                currentwidth += pls->size()+1;
            }
            out += *pls;
            out += space;
        }
    }

    //
    // display_descriptions returns a list of all of the class descriptions for a
    // Model and its child classes.
    //
    void
    display_descriptions(const Inheritance* inherit,const InheritanceList& ilist)
    {
        // Make sure the classes are registered...
        // Register_SCATMECH_Models();

        // Print out header...
        SCATMECH_output << endl << "-- Available " << inherit->get_name() << " Classes --" << endl;

        // Pointer to the list elements...
        InheritanceList::const_iterator v;

        // Point v to the begining of the list...

        int jj = 0,i;

        // Find the longest length of names...
        for (v = ilist.begin(); v != ilist.end(); ++v) {
            int k = ((*v)->get_name()).size();
            if ((*v)->children.size()>1 && (*v)!=inherit) {
                k+=7;
            }
            if (k>jj) jj=k;
        }

        for (i=0,v = ilist.begin(); v != ilist.end(); ++i,++v) {
            int k = ((*v)->get_name()).size();
            ostringstream message;
            message << "(" << (char)(i+'A') << ") ";
            SCATMECH_output << message.str();
            int menu=0;
            if ((*v)->children.size()>0 && (*v)!=inherit) {
                SCATMECH_output << "<MENU> ";
                menu=7;
            }
            SCATMECH_output << (*v)->get_name();
            if (Inheritance::show_model_mode) {
                for (int j=0; j<jj-k+2-menu; ++j) SCATMECH_output << '.';
                string message;
                Wrap_String((*v)->get_description(),message,jj+6,80);
                SCATMECH_output << message;
            }
            SCATMECH_output << endl;
        }
    }

    //
    // get_model displays a list of a Model's child classes and asks the user to
    // chose among them...
    //
    // If top==true, then one is at the top of the class heirarchy.
    //
    Model*
    Inheritance::
    get_model(bool top) const
    {
        // First...make sure all the SCATMECH models are registered...
        // Register_SCATMECH_Models();

        while (1) {

            InheritanceList ilist;

            get_model_inheritance_list(ilist,2);

            // If there is only one element in the database, return that
            // element...
            if (ilist.size()==1) return ilist[0]->make();

            InheritanceList::iterator v;
            for (v=ilist.begin(); v!=ilist.end(); ++v) {
                // Remove virtual elements that have no children or are the same as this...
                if ((*v)->is_virtual()&&((*v)->children.size()==0||(*v)==this)) {
                    v = ilist.erase(v);
                }
            }

            // Print out the class names...
            display_descriptions(this,ilist);

            // Print out some other options...
            SCATMECH_output << endl;
            if (show_model_mode) {
                SCATMECH_output << "(?) Turn off descriptions" << endl;
            } else {
                SCATMECH_output << "(?) Turn on descriptions" << endl;
            }

            if (!top) SCATMECH_output << "(^) Return to previous menu" << endl;

            string response = AskUser("\nChoose one",string("A"));

            Model* model = get_named_model(response,true);

            if (model) {
                return model;
            } else if (response[0]=='?') {
                show_model_mode= !show_model_mode;
            } else if (response[0]=='^' && !top) {
                return (Model*)0;
            } else if (response.size()==1) {
                // Convert alphanumeric response to index number...
                int k = toupper(response[0])-'A';
                // Make sure it is within range...
                if (k>=0&&k<(int)ilist.size()) {
                    // Get that model's instance...
                    InheritanceList::const_iterator v = ilist.begin();
                    advance(v,k);
                    const Inheritance* ib = *v;
                    if (ib==this) {
                        return ib->make();
                    } else {
                        Model* model = ib->get_model(false);
                        // If okay... were done...
                        if (model) {
                            return model;
                        } else {
                            return 0;
                        }
                    }
                }
            }
        }
    }

    //
    // Inheritance::get_parameters gets a list of a model's parameters...
    //
    void
    Inheritance::
    get_parameters(ModelParameterList& result,bool top) const
    {
        if (top) result.clear();
        if (parent) {
            parent->get_parameters(result,false);
        }
        for (ModelParameterList::const_iterator p=parameters.begin(); p!=parameters.end(); ++p) {
            result.push_back(*p);
        }
    }

    //
    // Inheritance::get_modelparameter searches the list of parameters associated with
    // the model, and returns the information about the parameter...
    //
    const  ModelParameterBase*
    Inheritance::get_modelparameter(const std::string& param) const
    {
        ModelParameterList mpl = get_parameters();
        for (ModelParameterList::const_iterator p=parameters.begin(); p!=parameters.end(); ++p) {
            if ((*p)->name == param) return (*p);
        }
        return NULL;
    }

    //
    // Inheritance::get_model_inheritance_list gets the list of all models inherited by
    // a specific model...
    //
    void
    Inheritance::
    get_model_inheritance_list(InheritanceList& result,int option) const
    {
        if (option) {
            result.clear();
            result.push_back(this);
        }

        for (InheritanceList::const_iterator v=children.begin(); v != children.end(); ++v) {
            result.push_back(*v);
            if (option!=2) (*v)->get_model_inheritance_list(result,0);
        }
    }

    //
    // Inheritance::get_named_model returns an instantiation of a class named model.
    //
    Model*
    Inheritance::
    get_named_model(const std::string& model,bool nothrow) const
    {
        const Inheritance *inheritance = get_named_inheritance(model,nothrow);
        if (inheritance) {
            if (!inheritance->virt)    return inheritance->make();
        }
        if (!nothrow) throw SCATMECH_exception("Model " + model + " is not instantiable");
        return (Model*)0;
    }

    const Inheritance*
    Inheritance::
    get_named_inheritance(const std::string& model,bool nothrow) const
    {
		if (model == name) return this;

        InheritanceList ilist;
        get_model_inheritance_list(ilist);

        for (InheritanceList::const_iterator v=ilist.begin(); v!=ilist.end(); ++v) {
            // If that model's name matches, return its inheritance...
            if (model==(*v)->get_name()) return *v;
        }

        // It was never found if it got here...
        if (!nothrow) {
            ostringstream os;
            os << string("Cannot find a " + name + " named " + model) << endl;
            os << "Valid choices are" << endl;
            for (InheritanceList::const_iterator v = ilist.begin(); v!=ilist.end(); ++v) {
                os << "   " << (*v)->get_name() << endl;
            }
            throw SCATMECH_exception(os.str());
        }  else {
            return (const Inheritance* )0;
        }
	} 

    ModelParameterBase::
    ModelParameterBase(const std::string& _name,
                       const std::string& _description,
                       const std::string& _type,
                       const std::string& _defaultvalue,
                       int _recalclevel,
                       Inheritance& inherit)
        : ParameterInfo(_name,_description,_type,_defaultvalue,_recalclevel)
    {
        inherit.add_parameter(this);
    }

    void
    Model::
    set_parameter_base(const std::string& parameter,const std::string& value)
    {
        if (parameter=="recalc") {
            recalc |= from_string<int>(value);
        } else if (parameter=="quiet") {
            int _quiet = from_string<int>(value);
            set_quiet(_quiet);
        } else {
            ModelParameterList mplist = get_inheritance().get_parameters();

            string _parameter = parameter;
            string _subparameter="";
            int pos = _parameter.find('.');

            if (pos!=string::npos) {
                _parameter = parameter.substr(0,pos);
                _subparameter = parameter.substr(pos+1,parameter.size());
                if (_subparameter=="")
                    error("Subparameter missing: " + parameter);
            }

            ModelParameterList::iterator q;
            for (q = mplist.begin(); q != mplist.end(); ++q) {
                const ModelParameterBase *param = *q;
                string name = param->name;

                if (_parameter==name) {
                    if (value!="default") {
                        param->set_parameter(this,_subparameter,value);
                    } else {
                        if (param->defaultvalue!="") {
                            param->set_parameter(this,_subparameter,param->defaultvalue);
                        }
                    }
                    return;
                }
            }

            error("Unknown parameter " + parameter + " in set_parameter_base");
        }
    }

    string
    Model::
    get_parameter_base(const string& parameter) const
    {
        ostringstream result;

        if (parameter=="recalc") {
            result << recalc;
        } else if (parameter=="quiet") {
            result << get_quiet();
        } else {
            ModelParameterList mplist = get_inheritance().get_parameters();

            string _parameter = parameter;
            string _subparameter="";
            int pos = _parameter.find('.');

            if (pos!=string::npos) {
                _parameter = parameter.substr(0,pos);
                _subparameter = parameter.substr(pos+1,parameter.size());
                if (_subparameter=="")
                    error("Subparameter missing: " + parameter);
            }

            ModelParameterList::iterator q;
            for (q = mplist.begin(); q != mplist.end(); ++q) {
                const ModelParameterBase *param = *q;
                string name = param->name;

                if (_parameter==name) {
                    result << param->get_parameter(this,_subparameter);
                    return result.str();
                }
            }

            error("Unknown parameter " + parameter + " in get_parameter");
        }
        return result.str();
    }

    ParameterInfo
    Model::
    get_parameter_info(const string& parameter) const
    {
        ParameterInfo result;
        if (parameter=="recalc") {
            result.name="recalc";
            result.type="int";
        } else if (parameter=="quiet") {
            result.name="quiet";
            result.type="int";
        } else {
            ModelParameterList mplist = get_inheritance().get_parameters();

            string _parameter = parameter;
            string _subparameter="";
            int pos = _parameter.find('.');

            if (pos!=string::npos) {
                _parameter = parameter.substr(0,pos);
                _subparameter = parameter.substr(pos+1,parameter.size());
                if (_subparameter=="")
                    error("Subparameter missing: " + parameter);
            }

            ModelParameterList::iterator q;
            for (q = mplist.begin(); q != mplist.end(); ++q) {
                const ModelParameterBase *param = *q;
                string name = param->name;

                if (_parameter==name) {
                    if (_subparameter.size()==0) {
                        result = *param;
                    } else {
                        if (param->get_inheritance(NULL)!=NULL) {
                            Model* m =(Model*)(param->get_ptr(this));
                            result = m->get_parameter_info(_subparameter);
                        } else {
                            error("Parameter " + parameter + " does not have subparameters");
                        }
                    }
                    return result;
                }
            }

            error("Unknown parameter " + parameter + " in get_parameter_info");
        }
        return result;
    }

    void
    Model::
    init()
    {
        ModelParameterList mplist = get_inheritance().get_parameters();

        ModelParameterList::iterator q;
        for (q = mplist.begin(); q != mplist.end(); ++q) {
            const ModelParameterBase *param = *q;
            if (param->defaultvalue!="") {
                param->set_parameter(this,"",param->defaultvalue);
            }
        }
    }

    void
    Model::
    AskUser()
    {
        init();

        ModelParameterList mplist = get_inheritance().get_parameters();

        ModelParameterList::iterator q;
        for (q = mplist.begin(); q != mplist.end(); ++q) {
            const ModelParameterBase* param = *q;
            param->AskUser(this);
        }
    }

    void
    Model::
    get_parameter_names_base(StringList& plist, const string& prefix) const
    {
        ModelParameterList mpl = get_inheritance().get_parameters();

        for (ModelParameterList::iterator q=mpl.begin(); q!=mpl.end(); ++q) {
            plist.push_back(prefix+(*q)->name);
            if ((*q)->get_inheritance(this)!=NULL) {
                Model* m = (Model*)((*q)->get_ptr(this));
                m->get_parameter_names_base(plist,prefix+(*q)->name+'.');
            }
        }
    }

    void
    Model::
    print_parameters_base(ostream& os,const string& prefix) const
    {
        //StringList vlist;
        //get_parameter_names(vlist);
        //StringList::iterator p;
        //for (p=vlist.begin();p!=vlist.end();++p) {
        //    ParameterInfo info = get_parameter_info(*p);
        //    os << *p << " = " << get_parameter(*p)
        //       << " (" << info.type << ": "
        //       << info.description << ')' << endl;
        //}
        ModelParameterList mplist = get_inheritance().get_parameters();
        for (ModelParameterList::iterator p = mplist.begin(); p!=mplist.end(); ++p) {
            const ModelParameterBase *param = *p;
            string name = param->name;
            string type = param->type;
            string description = param->description;
			string newprefix = (prefix != "") ? (prefix + ".") : ("");

            if (param->get_inheritance(this)==NULL) {
                os << newprefix << name << " = " << param->get_parameter(this,"")
                   << " (" << type << ": "
                   << description << ")" << endl;
            } else {
                os << newprefix << name << " = " << param->get_parameter(this,"")
                   << " (" << type << ": "
                   << description << ")" << endl;

                ((Model*)(param->get_ptr(this)))->print_parameters_base(os,newprefix+name);
            }
        }
    }

    void
    Model::
    error(const string& message) const
    {
        string message2 = "Error in class " + get_inheritance().get_name() + ": " + message + '\n';
        ostringstream out;
        print_parameters(out);
        throw SCATMECH_exception(message2+out.str());
    }

    Inheritance Model::inheritance("Model","Generalized Model","",0,0);


} // namespace SCATMECH


