//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: inherit.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_INHERIT_H
#define SCATMECH_INHERIT_H

#include <deque>
#include "scatmech.h"
#include <memory>

#define SMERROR(message)    error(message,__FILE__,__LINE__);

namespace SCATMECH {

    //
    // This file defines, amongst other things, the following classes...
    //
    class Inheritance;
    class Model;
    template <class _Model> class Model_Ptr;
    template <class _Model> class Model_Maker;
    class Model_Maker_Base;
    struct ParameterInfo;
    class ModelParameterBase;
    template <class MODEL,class TYPE> class ModelParameter;
    template <class MODEL,class TYPE> class ModelPtrParameter;
    class Get_Model_Ptr;

    template <class _Model> _Model* _Get_Named_Model(const STRING& name);
    template <class _Model> _Model* _Get_Model();

    typedef std::deque<const Inheritance*> InheritanceList;
    typedef std::unique_ptr<ModelParameterBase> ModelParameterBasePtr;
    typedef std::deque<const ModelParameterBase*> ModelParameterList;
    typedef std::deque<STRING> StringList;

    ///
    /// The struct ParameterInfo stores information about a parameter.
    ///
    struct ParameterInfo {
        /// Default constructor.
        ParameterInfo() {}
        /// Constructor
        ParameterInfo(
            const STRING& _name,         ///< The name of the parameter
            const STRING& _description,  ///< Its description
            const STRING& _type,         ///< The data type
            const STRING& _defaultvalue, ///< Its default value
            int _recalclevel             ///< Its recalc level
        ) : name(_name), description(_description),
            type(_type), defaultvalue(_defaultvalue),
            recalclevel(_recalclevel) {}

        STRING name;           ///< The name of the parameter...
        STRING description;    ///< A description of the parameter...
        STRING type;           ///< The data type for external programs...
        STRING defaultvalue;   ///< The default value for the parameter...
        int recalclevel;       ///< A flag defining what needs to be done
        ///< when this parameter has been modified (usually -1)...
    };

    ///
    /// @brief Base class for all Model parameters.
    ///
    /// It stores the pointer to a member in a Model.
    ///
    class ModelParameterBase : public ParameterInfo {
        public:
            /// Constructor
            ModelParameterBase(
                const STRING& _name,         ///< The name of the parameter
                const STRING& _description,  ///< Its description
                const STRING& _type,         ///< The data type
                const STRING& _defaultvalue, ///< Its default value
                int _recalclevel,            ///< The recalc level
                Inheritance& inherit         ///< The inheritance structure to which
                ///< model this parameter belongs
            );

            /// Set a parameter's value.
            virtual void set_parameter(
                Model* model,               ///< Pointer to the specific object
                const STRING& subparameter, ///< Any subparameters, if applicable
                const STRING& value         ///< A string representation of the value
            ) const = 0;

            /// Return a pointer to the instance of the parameter.
            /// How this pointer is interpreted should depend upon type.
            virtual void* get_ptr(
                const Model* model       ///< Pointer to the specific object
            ) const = 0;

            /// Return the current value of the parameter as a STRING.
            virtual STRING get_parameter(
                const Model* model,         ///< Pointer to the specific object
                const STRING& subparameter  ///< Any subparameters, if applicable
            ) const=0;

            /// Prompts the user for the parameter's value.
            virtual void AskUser(
                Model* model        ///< Pointer to the specific object
            ) const = 0;

            /// Return parent class' inheritance.  Function returns a null pointer
            /// for a normal parameter, and a pointer the parent class'
            /// inheritance if it is a Model_Ptr parameter...
            virtual const Inheritance* get_inheritance(
                const Model* model=NULL ///< Pointer to the specific object
            ) const {
                return (Inheritance*)0;
            }
    };

    ///
    /// @brief Base class for all models in SCATMECH
    ///
    class Model
    {
        public:
            /// Constructor
            Model() {
                recalc=-1;
                init();
            }

            /// Destructor
            virtual ~Model() {}

            /// @brief Initialize the model.
            /// Assigns all the models parameters to their
            /// default values.
            void init();

            /// @brief Query the user for parameters.
            /// A default implementation if this function is
            /// provided, which can be used under most conditions.
            void AskUser();

            /// Static variable that keeps track of all inherited models.
            static Inheritance inheritance;

            /// @brief Return the inheritance of the specific model.
            /// It should be redefined for every inherited model, so that it will
            /// return the appropriate model.  The DECLARE_MODEL macro does this.
            virtual const Inheritance& get_inheritance() const {
                return inheritance;
            }

            /// Return a list of all the parameters.
            void get_parameter_names(
                StringList& plist ///< Where the list of parameters is returned
            ) const {
                get_parameter_names_base(plist,"");
            }

            /// Return information about a specific parameter.
            ParameterInfo get_parameter_info(
                const STRING& parameter ///< The name of the parameter
            ) const;

            /// Set a parameter to a value.
            template <class SETTYPE>
            void set_parameter(
                const STRING& parameter, ///< The name of the parameter
                const SETTYPE& value     ///< The value to set it to
            ) {
                set_parameter_base(parameter,to_string(value,16));
            }

            /// Set a parameter to a value.
            void set_parameter(
                const STRING& parameter, ///< The name of the parameter
                const char* value        ///< The value to set it to
            ) {
                set_parameter_base(parameter,std::string(value));
            }

            /// Get the current value of a parameter in string form
            STRING get_parameter(
                const STRING& parameter ///< The parameter name
            ) const {
                return get_parameter_base(parameter);
            }

            /// @brief Output parameters to a stream.
            /// Send a listing of all parameters and their current values to an output stream.
            void print_parameters(
                std::ostream& os,               ///< Output stream
                const std::string& prefix = ""  /// Prefix (not including ".") to be placed before each parameter name
            ) const {
                print_parameters_base(os,prefix);
            }

            /// Set by OR'ing to force recalculation.
            void set_recalc(int _recalc=0xFF) {
                recalc |= _recalc;
            }

            /// Returns non-zero if any parameter has changed.  This member function should
            /// be overloaded if other changes require setup to be called.
            virtual int get_recalc() const {
                return recalc;
            }

            /// Set to keep model quiet during calculations.
            static void set_quiet(int _quiet) {
                static mutex_t mutex;
                mutex.lock();
                quiet=_quiet;
                mutex.unlock();
            }

            /// Returns 1 if model should keep quiet.
            static int get_quiet() {
                return quiet;
            }

            /// @brief Perform housekeeping, if neccessary.
            /// All routines that need up-to-date information should begin with SETUP().
            void SETUP() {
                if (get_recalc()) setup();
            }

        protected:
            /// @brief Routine which does housekeeping if recalc!=0.
            /// All inherited classes which need to perform calculations
            /// once anytime a parameter changes should define setup(),
            /// and should begin by calling its parent's setup().
            virtual void setup() {
                RECALC = recalc;
                recalc=0;
            }

            /// Routine to send messages to the user,
            void message(
                const STRING& s  ///< Message to be sent
            ) const {
                if (!quiet) {
                    SCATMECH_output << s << ret;
                }
            }

            /// Routine which handles errors...
            void error(const STRING& message) const;

            /// @brief Protected virtual version of set_parameter().
            /// This parameter is a non-templated version which handles
            /// the core functions of set_parameter. Any class defining this function
            /// should call its parent's set_parameter_base function, if it cannot find the parameter.
            virtual void set_parameter_base(
                const STRING& parameter, ///< The parameter name
                const STRING& value      ///< String represention of a value
            );

            /// @brief Gets a parameter to a value.
            /// This parameter is a non-templated version which handles
            /// the core functions of set_parameter.  Any class defining this function
            /// should call its parent's set_parameter_base function, if it cannot find the parameter.
            virtual STRING get_parameter_base(
                const STRING& parameter ///< The parameter name
            ) const;

            /// @brief get_parameter_names(plist)
            /// This parameter is a non-templated version which handles
            /// the core functions of get_parameter_names. The parameter prefix
            /// allows for nested model names. Any class defining this function should
            /// call its parent's get_parameter_names_base before adding specialized parameters.
            virtual void get_parameter_names_base(
                StringList& plist,      ///< Receives result of query
                const STRING& prefix    ///< Prefix to add to resulting parameter names
            ) const;

            /// @brief Output parameters to a stream.
            /// Virtual version of print_parameters.  Send a listing of all parameters and their current values to an output stream.
            /// Any class defining this function should call its parent's print_parameter_names_base
            /// before adding specialized parameters.
            virtual void print_parameters_base(
                std::ostream& os,               ///< Output stream
                const std::string& prefix = ""  /// Prefix (not including ".") to be placed before each parameter name
            ) const;


            /// @brief During setup(), RECALC holds the current value of recalc.
            /// This protected variable can be used to determine which steps need to
            /// be taken during setup().
            int RECALC;

        private:

            /// @brief Indicates that initialization may need to be carried out.
            /// Should be set anytime a variable is set.
            int recalc;

            /// @brief Indicates model should be verbose.
            /// Indicates whether or not the model should be verbose during its
            /// calculations.  By default, quiet is zero.
            static int quiet;

    };

    /// @brief Class that keeps track of class inheritance structure.
    ///
    /// class Inheritance keeps track of all of the child classes of a given
    /// parent class.   It enables a specific instance of a child class to be
    /// created by a menu or by a character string.  It keeps track of a name,
    /// description, parametres, and where the Model is in the Model heirarchy.
    /// It provides functions for linking to its parents and children and
    /// progeny (children's children).
    ///
    class Inheritance {
        public:
            /// Constructor...
            Inheritance(
                const STRING& _name,      ///< Name of class
                const STRING& _desc,      ///< Description
				const STRING& _parentname, ///< Parent's name
                Model_Maker_Base *_maker, ///< Class which creates new instances of class
                Inheritance *_parent      ///< Pointer to Parent inheritance
            ) : name(_name),description(_desc),parentname(_parentname),
                maker(_maker),parent(_parent),
                virt(_maker==0),registered(false) {}

            /// Return name of the class.
            const STRING& get_name() const {
                return name;
            }

            /// Return desscription of the class.
            const STRING& get_description() const {
                return description;
            }

            /// Return true if it is a virtual model
            bool is_virtual() const {
                return virt;
            }

            /// Return true if model can be instantiated.
            bool is_instantiable() const {
                return !virt;
            }

            /// @brief Query for child class.
            /// Print out a list of available child classes and
            /// return a pointer to the user's choice.
            Model* get_model() const {
                return get_model(true);
            }

            /// Return a pointer to the named progeny class.
            Model* get_named_model(
                const STRING& name   ///< Name of class to create instance of
            ) const {
                return get_named_model(name,false);
            }

            /// Returns parent's Inheritance.
            const Inheritance* get_parent() const {
                return parent;
            }

            /// @brief Return list of direct children.
            /// Direct children includes only directly inherited children.
            const InheritanceList& get_children() const {
                return children;
            }

            /// @brief Return list of all children.
            /// Progeny includes children, children's children, etc.
            InheritanceList get_progeny() const {
                InheritanceList result;
                get_model_inheritance_list(result,1);
                return result;
            }

            /// Return the Inheritance for a named progeny.
			const Inheritance*
			get_named_inheritance(
				const STRING& model,  ///< Name of the model
				bool nothrow = false    ///< Will throw exception if not found.
			) const;

            /// Get list of model parameters.
            ModelParameterList get_parameters() const {
                ModelParameterList result;
                get_parameters(result,true);
                return result;
            }

            /// Get a specific model parameter.
            const ModelParameterBase* get_modelparameter(
                const STRING& param   ///< Name of parameter
            ) const;

            // Return an instance of the model.
            Model* make() const;

            /// Return a clone of the model.
            Model* clone(
                const Model& m   ///< Model to clone
            ) const;

            /// @brief Register the class' existence.
            /// Registers with its own list and with its parent's list.
            void Register_Model(Inheritance& topModel = Model::inheritance) {
                if (!registered) {
					if (parent) { 
						// This part was modified 6/2018 in a step to enable models to be added by DLL. 
						// parent->add_child(this);
						Inheritance* in = const_cast<Inheritance*>(topModel.get_named_inheritance(parentname));
						if (in) in->add_child(this);
					}
                    registered=true;
                }
            }


        private:
            friend class ModelParameterBase;

            /// Adds a class to the list this class' inherited classes...
            void add_child(
                Inheritance* ib   ///< Pointer to class inheritance to add
            ) {
                children.push_back(ib);
            }

            /// Add a model parameter to this model.
            void add_parameter(
                ModelParameterBase* mp   ///< Pointer to parameter
            ) {
                parameters.push_back(mp);
            }

            /// A string name for the class
            STRING name;

            /// A string description for the class
            STRING description;

			/// A string giving the name of the class' parent
			STRING parentname;

            /// A list of pointers to inherited classes
            InheritanceList children;

            /// Points to the Model's parent's Inheritance.
            Inheritance* parent;

            /// The list of parameters for this model.
            ModelParameterList parameters;

            /// A flag to indicate if this class is instantiable or not...
            bool virt;

            /// A flag to indicate if this class has been registered.
            bool registered;

            /// A static flag that indicates whether or not the user
            /// would like to see descriptions with the class names.
            static bool show_model_mode;

            /// The following contains a pointers to the class that contain the
            /// maker and the registrars for this Model class...
            std::unique_ptr<Model_Maker_Base> maker;

            //
            // Functions that returns a list of inherited classes and their
            // descriptions....
            //
            // option = 1 indicates result should be cleared at start
            // option = 2 indicates that only first level should be traversed

            void get_model_inheritance_list(InheritanceList& result,int option=1) const;
            Model* get_model(bool top) const;
            Model* get_named_model(const STRING& name,bool nothrow) const;
			void get_parameters(ModelParameterList& result, bool top) const;
			friend void display_descriptions(const Inheritance* inherit,
                                             const InheritanceList& ilist);
    };

    /// @brief Smart pointer to a Model
    ///
    /// The template class Model_Ptr<_Model> works similarly to std::unique_ptr<T>,
    /// except that it makes cloned copies rather than transferring ownership of
    /// a pointer.  It also adds some functions which make use of the SCATMECH
    /// Inheritance and Model classes.
    template <class _Model>
    class Model_Ptr
    {
        public:
            typedef _Model Type;
            typedef _Model* Ptr;

            /// Default constructor creates an uninitialized pointer.
            Model_Ptr() : model(0) {}

            /// Copy constructor
            Model_Ptr(const Model_Ptr<_Model>& model_ptr) {
                if (model_ptr.model) {
                    const Inheritance& i = model_ptr.model->get_inheritance();
                    model = static_cast<_Model*>(i.clone(*(model_ptr.model)));
                }
                else  model = 0;
            }

            /// Assignment operator
            Model_Ptr& operator=(const Model_Ptr<_Model>& model_ptr) {
                if (model) delete model;
                if (model_ptr.model) {
                    const Inheritance& i = model_ptr.model->get_inheritance();
                    model = static_cast<_Model*>(i.clone(*(model_ptr.model)));
                }
                else model=0;
                return *this;
            }

            /// @brief Constructor with a pointer
            /// Ownership is transferred.
            Model_Ptr(Model* m) : model(0) {
                *this = m;
            }

            /// @brief Assignment to a pointer
            /// Ownership is transferred.
            Model_Ptr& operator=(Model* m) {
                if (model) delete model;
                if (m) {
                    const Inheritance& i = m->get_inheritance();
                    model = static_cast<_Model*>(m);
                }
                else model = 0;
                return *this;
            }

            /// @brief Query the user for an inherited model.
            /// The constructor with the trivial class Get_Model_Ptr
            /// gets an instance of the class using the
            /// Inheritance functionality.
            Model_Ptr(const Get_Model_Ptr& gm) : model(0) {
                *this = gm;
            }

            /// @brief Query the user for an inherited model.
            /// The assignment to the trivial class Get_Model_Ptr
            /// gets an instance of the class using the
            /// Inheritance functionality.
            Model_Ptr& operator=(const Get_Model_Ptr& gm) {
                GetPtr();
                return *this;
            }

            /// @brief Assignment to a string.
            /// Gets an instance of the class using the Inheritance functionality which
            /// has a specific name.
            Model_Ptr& operator=(const STRING& name) {
                if (model) delete model;
                model = _Get_Named_Model<_Model>(name);
                return *this;
            }

            /// @brief Constructor with a string.
            /// Gets an instance of the class using the Inheritance functionality which
            /// has a specific name.
            Model_Ptr(const STRING& name) : model(0) {
                *this = name;
            }

            /// @brief Destructor
            /// Destructor deletes the pointer if it is valid.
            ~Model_Ptr() {
                if (model) delete model;
            }

            /// Access elements of model as if Model_Ptr<_Model> were a normal pointer...
            _Model* operator->() const {
                if (!model) throw(SCATMECH_exception(
                                          "Attempt to access element by null pointer"));
                return get();
            }

            /// Dereferencing
            _Model& operator*() const {
                return *(get());
            }

            /// Return the pointer
            _Model* get() const {
                return model;
            }

            /// Return an instance of _Model using _Model's Inheritance structure...
            Model_Ptr& GetPtr() {
                if (model) delete model;
                model =  _Get_Model<_Model>();
                return *this;
            }

            /// Returns the model's inheritance.
            static const Inheritance* GetInheritance() {
                return &(_Model::inheritance);
            }

        private:
            /// The pointer to a Model...
            _Model* model;
    };


    /// @brief Trivial class used by Model_Ptr<_Model>
    ///
    /// Trivial class used by the constructor and operator=() of
    /// Model_Ptr<_Model> to indicate that get_model must be called...
    class Get_Model_Ptr
    {};

    /// @brief Class for creating and cloning models
    ///
    /// The templated class Model_Maker<_Model> is used by the class Inheritance
    /// to provide a function which creates a new instance of a class or a clone
    /// of a class.
    class Model_Maker_Base {
        public:
            /// @brief Create an instance of the model
            /// It is expected that the programmer will typecast the
            /// result to a pointer to the approprate model.
            virtual Model* make() const = 0;

            /// @brief Create a clone of the model
            /// It is expected that the programmer will typecast the
            /// result to a pointer to the approprate model.
            virtual Model* clone(const Model& model) const = 0;
    };

    ///
    /// @brief Template class to create and clone specific models
    ///
    template <class _Model>
    class Model_Maker : public Model_Maker_Base {
        public:
            Model* make() const {
                Model* newModel = static_cast<Model*>(new _Model);
                newModel->init();
                return newModel;
            }
            Model* clone(const Model& model) const {
                return new _Model(static_cast<const _Model&>(model));
            }
    };

    ///
    /// @brief Placeholder template for something which is not a model
    ///
    template <>
    class Model_Maker<void> : public Model_Maker_Base {
        public:
            Model* make() const {
                return static_cast<Model*>(0);
            }
            Model* clone(const Model& model) const {
                return static_cast<Model*>(0);
            }
    };

    inline Model* Inheritance::make() const {
        return maker->make();
    }
    inline Model* Inheritance::clone(const Model& m) const {
        return maker->clone(m);
    }

    /// @brief Register all models in the tree of classes
    ///
    /// Model that have children define this function, so that the tree
    /// can be built.  This function should be called before any
    /// calls that search the class tree.
    void Register(
        const Model* x    ///< Argument used for defining polymorphic function
    );

    /// @brief Returns a pointer to a _Model.
    ///
    /// Queries the user for a model, creates an instance of one, and returns
    /// the pointer to it.  It registers the model first.
    template <class _Model>
    _Model* _Get_Model()
    {
        Register((_Model*)0);
        return static_cast<_Model*>(_Model::inheritance.get_model());
    }

    /// @brief Returns a pointer to a specified Model.
    ///
    /// Creates an instance of the model with the given name.
    /// It registers the model first.
    template <class _Model>
    _Model* _Get_Named_Model(
        const STRING& name  ///< The name of the model (must be a child of _Model)
    ) {
        Register((_Model*)0);
        return static_cast<_Model*>(_Model::inheritance.get_named_model(name));
    }

    /// Macro to register a specific model
#define Register_Model(_Model) _Model::inheritance.Register_Model()

    //#define Get_Model(_Model)     SCATMECH::_Get_Model<_Model>()
    //#define Get_Named_Model(_Model,name)  SCATMECH::_Get_Named_Model<_Model>(name)

    /// Routine to assert that a string is empty
    inline void
    AssertNullParameter(const STRING& parameter)
    {
        if (parameter!="") throw SCATMECH_exception("No subparameter allowed: " + parameter);
    }

    /// @brief Sets a parameter to a value, given a string.
    ///
    /// Used by ModelParameter, the template function ModelParameterSet converts a string
    /// to a TYPE, using subparameter if the data type handles subparameters. By default,
    /// it does not allow for non-empty subparameters.  New data types, which act as
    /// model parameters, need to define an explicit template specialization, if it
    /// allows for subparameters.
    template <class TYPE>
    void ModelParameterSet(
        TYPE& variable,                 ///< A reference to the variable being set
        const STRING& subparameter,     ///< The subparameter to set (default asserts empty)
        const STRING& value             ///< String representation of new value
    )
    {
        AssertNullParameter(subparameter);
        variable = from_string<TYPE>(value);
    }

    /// @brief Gets a parameter specified by a string.
    ///
    /// Used by ModelParameter, the template function ModelParameterGet converts TYPE to
    /// STRING, using subparameter if the data type handles subparameters. By default,
    /// it does not allow for non-empty subparameters.  New data types, which act as
    /// model parameters, need to define an explicit template specialization, if it
    /// allows for subparameters.  It returns a string representation of the parameter.
    template <class TYPE>
    STRING ModelParameterGet(
        TYPE& variable,                 ///< Reference to the variable
        const STRING& subparameter      ///< The subparameter to set (default asserts empty)
    )
    {
        AssertNullParameter(subparameter);
        return to_string(variable);
    }

    /// @brief Query the user for a parameter value
    ///
    /// Used by ModelParameter, the template function ModelParameterAskUser sends
    /// a query to the user and sets the value to the response. New data types,
    /// which act as model parameters, need to define an explicit template
    /// specialization, if it has subparameters.
    template <class TYPE>
    void ModelParameterAskUser(
        TYPE& variable,         ///< Reference to the variable
        const STRING& prompt    ///< The prompt used to query
    )
    {
        variable = AskUser(prompt,variable);
    }

    /// @brief Simple (non-Model_Ptr) parameters
    ///
    /// Simple (non-Model_Ptr) parameters use this class
    template <class MODEL,class TYPE>
    class ModelParameter : public ModelParameterBase {
        public:
            ModelParameter(const STRING& _name,
                           const STRING& _description,
                           const STRING& _type,
                           const STRING& _defaultvalue,
                           int _recalclevel,
                           Inheritance& inherit,
                           TYPE MODEL::*_parameter)
                : ModelParameterBase(_name,_description,_type,_defaultvalue,_recalclevel,inherit),
                  parameter(_parameter) {}

            virtual void set_parameter(Model* model,
                                       const STRING& subparameter,
                                       const STRING& value) const
            {
                TYPE &_parameter = ((MODEL*)model)->*parameter;
                ModelParameterSet<TYPE>(_parameter,subparameter,value);
                model->set_recalc(recalclevel);
            }

            virtual void AskUser(Model* model) const {
                TYPE &_parameter = ((MODEL*)model)->*parameter;
                ModelParameterAskUser<TYPE>(_parameter,description);
                model->set_recalc(recalclevel);
            }

            virtual void* get_ptr(const Model* model) const {
                return (void*)&(((MODEL*)model)->*parameter);
            }

            virtual STRING get_parameter(const Model* model,
                                         const STRING& subparameter) const
            {
                TYPE &_parameter = ((MODEL*)model)->*parameter;
                return ModelParameterGet<TYPE>(_parameter,subparameter);
            }

        private:
            // Pointer to member variable...
            TYPE MODEL::*parameter;
    };

    /// @brief Model_Ptr parameters
    /// Model_Ptr parameters use this class
    template <class MODEL,class TYPE>
    class ModelPtrParameter : public ModelParameterBase {
        public:
            ModelPtrParameter(const STRING& _name,
                              const STRING& _description,
                              const STRING& _type,
                              const STRING& _defaultvalue,
                              int _recalclevel,
                              Inheritance& inherit,
                              TYPE MODEL::*_parameter)
                : ModelParameterBase(_name,_description,_type,_defaultvalue,_recalclevel,inherit),
                  parameter(_parameter) {}

            virtual void set_parameter(Model* model,
                                       const STRING& subparameter,
                                       const STRING& value) const
            {
                TYPE &_parameter = ((MODEL*)model)->*parameter;

                if (subparameter=="") {
                    // If there is no subparameter, then create an
                    // instance of the model type...
                    Register(TYPE().get());
                    _parameter = TYPE::GetInheritance()->get_named_model(value);
                } else {
                    // If there is a subparameter, then send it to the model...
                    _parameter->set_parameter(subparameter,value);
                }
                model->set_recalc(recalclevel);
            }

            virtual void AskUser(Model* model) const {
                TYPE &_parameter = ((MODEL*)model)->*parameter;
                SCATMECH_output << description;
                _parameter = Get_Model_Ptr();
				_parameter->Model::AskUser();
                model->set_recalc(recalclevel);
            }

            virtual const Inheritance* get_inheritance(const Model* model) const {
                if (model!=NULL) {
                    return &((((MODEL*)model)->*parameter).get()->get_inheritance());
                } else {
                    return TYPE::GetInheritance();
                }
            }

            virtual void* get_ptr(const Model* model) const {
                return (void*)((((MODEL*)model)->*parameter).get());
            }

            virtual STRING get_parameter(const Model* model,
                                         const STRING& subparameter) const
            {
                if (subparameter=="") {
                    return ((MODEL*)model->*parameter).get()->get_inheritance().get_name();
                } else {
                    return ((MODEL*)model->*parameter).get()->get_parameter(subparameter);
                }
            }

        private:

            // A pointer to member for the actual parameter
            TYPE MODEL::*parameter;
	};

	///
	/// @brief Macro used in the declaration of a Model
	///
	/// Each model should have this macro in its declaration.  It declares a
	/// static Inheritance for the model, a get_inheritance() function, and
	/// a clone() function.
	///
#define DECLARE_MODEL() \
		public: static SCATMECH::Inheritance inheritance; \
		public: SCATMECH::Model* clone() const {return inheritance.clone(*this);} \
		public: virtual const SCATMECH::Inheritance& get_inheritance() const \
						{return inheritance;}

	///
	///  @brief Macro declare a parameter in a class declaration.
	///
	///  Each parameter in a model should be declared using this macro.  It declares
	///  the parameter as protected, member set_ and get_ functions as public,
	///  and a static private ModelParameter variable.
	///
	///  @param TYPE is the data type for the parameter
	///  @param PARAMETER is the variable name.
	///
#define DECLARE_PARAMETER(TYPE,PARAMETER) \
		public:    const TYPE& get_##PARAMETER() const {return PARAMETER;} \
		public:    void set_##PARAMETER(const TYPE& __x) {PARAMETER = __x; set_recalc(MP_##PARAMETER->recalclevel);} \
		protected: TYPE PARAMETER; \
		private: static SCATMECH::ModelParameterBasePtr MP_##PARAMETER;

	///
	/// @brief Macro to define a instantiable model's static Inheritance element
	///
	/// DEFINE_MODEL aids the programmer in defining a Model's
	/// static Inheritance element.
	///
	/// @param child is the class which has the member inheritance which needs defining.
	/// @param parent is the parent class
	/// @param name is a string version of the name of the class
	/// @param desc is a string description of the class
	///
#define DEFINE_MODEL(child,parent,desc) \
	SCATMECH::Inheritance child::inheritance(#child, \
		desc, #parent, \
		new SCATMECH::Model_Maker<child>(), \
		&parent::inheritance);

    ///
    /// @brief Macro to define a non-instantiable model's static Inheritance element
    ///
    /// DEFINE_VIRTUAL_MODEL aids the programmer in defining a Model's
    /// static Inheritance element.
    ///
    /// @param child is the class which has the member inheritance which needs defining.
    /// @param parent is the parent class
    /// @param name is a string version of the name of the class
    /// @param desc is a string description of the class
    ///
#define DEFINE_VIRTUAL_MODEL(child,parent,desc) \
	SCATMECH::Inheritance child::inheritance(#child, \
		desc, #parent, \
		0, \
		&parent::inheritance);

    ///
    /// @brief Macro to define non-Model_Ptr parameters
    ///
    /// DEFINE_PARAMETER defines the static ModelParameter variable associated
    /// with a specific model parameter.
    ///
    /// This macros must follow the DEFINE_MODEL macro for CLASS.
    ///
    /// @param CLASS is the Model class
    /// @param TYPE is the data type for the parameter
    /// @param PARAMETER is the parameter name.
    /// @param DESCRIPTION is a string description of the parameter
    /// @param DEFAULTVALUE is a string representation of the default value
    ///
#define DEFINE_PARAMETER(CLASS,TYPE,PARAMETER,DESCRIPTION,DEFAULTVALUE,RECALCLEVEL) \
		ModelParameterBasePtr \
		CLASS::MP_##PARAMETER( \
							  new ModelParameter< CLASS , TYPE >(#PARAMETER, \
																 DESCRIPTION, \
																 #TYPE, \
																 DEFAULTVALUE, \
																 RECALCLEVEL, \
																 CLASS::inheritance, \
																 &CLASS::PARAMETER \
																) \
							  );

    ///
    /// @brief Macro to define Model_Ptr parameters
    ///
    /// DEFINE_PARAMETER defines the static ModelPtrParameter variable associated
    /// with a specific model Model_Ptr parameter.
    ///
    /// This macros must follow the DEFINE_MODEL macro for CLASS.
    ///
    /// @param CLASS is the Model class
    /// @param TYPE is the data type for the parameter
    /// @param PARAMETER is the parameter name.
    /// @param DESCRIPTION is a string description of the parameter
    /// @param DEFAULTVALUE is a string representation of the default value
    ///
#define DEFINE_PTRPARAMETER(CLASS,TYPE,PARAMETER,DESCRIPTION,DEFAULTVALUE,RECALCLEVEL) \
		ModelParameterBasePtr \
		CLASS::MP_##PARAMETER( \
							  new ModelPtrParameter< CLASS , TYPE >(#PARAMETER, \
																	DESCRIPTION, \
																	#TYPE, \
																	DEFAULTVALUE, \
																	RECALCLEVEL, \
																	CLASS::inheritance, \
																	&CLASS::PARAMETER \
																   ) \
							 );

} // namespace SCATMECH

#endif
