import SCATPY

class Model:
    """
    Base class for Models.
    """

    def __init__(self,*args,**kwargs):
        self.newer = SCATPY.Get_Model
        self.freer = SCATPY.Free_Model

        if len(args)>0 and type(args[0])==str:
                self.handle = self.newer(args[0])
                args = args[1::]
        else:
            self.handle = self.newer("Microroughness_BRDF_Model")
        self.setParameters(*args, **kwargs)

    def setParameters(self,*args, **kwargs):
        """
        Sets one or more parameters for the model. 

        Parameters
        ----------
        parameter_values : dict 
                           Argument is a dictionary with with str keys being 
                           the parameters and values being the values. 
                           `None` key defines model name.

        **kwargs : 
                    parameters can be set by parameter name, 
                    `lambda` should be replace with `wavelength`, since
                    lambda is a Python keyword.

        Returns
        -------
        Model


        """
        # Change `wavelength` to `lambda` in kwargs...
        if "wavelength" in kwargs:
            kwargs["lambda"]=kwargs["wavelength"]
            del kwargs["wavelength"]

        parameter_values = dict()
        for arg in args:
            if hasattr(arg,"getParameters"):
                arg = arg.getParameters()
            parameter_values.update(arg)
            
        parameter_values.update(kwargs)
            
        for parameter,value in parameter_values.items():
            
            if hasattr(value,"getModelDict"):
                # value is a FilmStack...
                value = value.getModelDict()
            elif hasattr(value,"nkstr"):
                # value is a OpticalFunction...
                value = value.nkstr()
            elif type(value) is complex:
                # value is complex and needs to be converted to a str...
                value = "(%16g,%16g)" % (value.real,value.imag)

            if parameter==None:
                # Set the Model...
                self.freer(self.handle)
                self.handle = self.newer(value)
            elif type(value)==dict:
                dictnew = {}
                for subparameter,subvalue in value.items():
                    if hasattr(subvalue,'imag'):
                        if subvalue.imag==0:
                            subvalue = "%.16g" %subvalue.real
                        else:
                            subvalue = "(%.16g,%.16g)" % (subvalue.real,subvalue.imag)
                        
                    if type(subparameter)==type(None):
                        SCATPY.SetParameter(self.handle,
                                            str(parameter),
                                            str(subvalue))
                    else:
                        dictnew[parameter + "." + subparameter] = subvalue
                self.setParameters(dictnew)
            elif isinstance(value,Model):
                self.setParameters({str(parameter) : value.getParameters()})
            else:
                SCATPY.SetParameter(self.handle,str(parameter),str(value))
        return self

    def setParameter(self, parameter, value):
        """
        Sets a single parameter for the model.

        Arguments
        ---------
        parameter: str
                   Parameter name as a string.  
                   Subparameters must use `parameter.subparameter` notation.
        value: float, str, or complex
               Parameter value
        """
        return self.setParameters({parameter : value})
        
        if hasattr(value,'imag'):
            if value.imag==0:
                value = "%.16g" % value.real
            else:
                value = "(%.16g,%.16g)" % (value.real,value.imag)

        SCATPY.SetParameter(self.handle,str(parameter),value)
        return self
            
    def getParameter(self, parameter):
        """
        Gets a string corresponding to the value of a model parameter.
        Arguments:
        parameter: parameter name as a string
        Returns: parameter value
        """
        return SCATPY.GetParameter(self.handle, str(parameter))

    def getParameters(self):
        """
        Gets all of the parameters for the model.

        Returns
        -------
             dict with str keys
        """
        dict = SCATPY.GetParameterDictionary(self.handle)
        newdict = {None : self.getModelName()}
        for d in dict:
            newdict[d]=dict[d]['value']
        return newdict
        
    def printParameters(self):
        """
        Returns a printable string that includes all the parameters, 
        their values, and their descriptions.
        Recommend surround with print()
        """
        return SCATPY.PrintParameters(self.handle)

    def getParameterDictionary(self):
        """
        Returns a dictionary of dictionaries for each parameter that includes 
        values, descriptions, and types.
        """
        return SCATPY.GetParameterDictionary(self.handle)

    def askParameters(self):
        """
        Queries the user for parameter values.
        TODO FIX!
        """ 
        params = self.getParameterDictionary()
        for param in params.copy():
            if '.' in param:
                del params[param]
        del param
        result = {}
        for param,d in params.items():
            if d['type'][-4:]!="_Ptr":
                response = input(d['description']+" <"+ d['value'] + ">:")
                if response == '':
                    result[param] = str(d['value'])
                else:
                    result[param] = response
            else:
                result[param]= chooseModelAndParameters(d['type'][:-4])
        self.setParameters(result)

    def getModelName(self):
        """
        Returns the model's name
        """
        return SCATPY.GetModelName(self.handle)

    def __repr__(self):
        return self.getModelName()
    
    def __str__(self):
        """
        Returns a string version of the model, which is a pretty-printed map 
        of maps.
        """
        result = {None : self.getModelName()}
        params = self.getParameters()
        for p in params:
            result[p] = params[p]
        return strDictOfDict(createDictOfDict(result))
        
    def __del__(self):
        try:
            SCATPY.Free_Model(self.handle)
        except:
            pass
                
def getModelDictionary(model="Model"):
    """
    Returns a dictionary containing all the models inheriting model and their 
    descriptions.
    """
    return SCATPY.GetModelDictionary(model)

def chooseModel(model):
    """
    Presents user with a choice of models and returns the user's choice, by 
    name.
    """
    choices = getModelDictionary(model)
    i = 0
    for m in choices:
        choices[m][None] = i
        i=i+1
        
    valid = False
    while not valid:
        i = 0
        for m,d in choices.items():
            print(d[None]," : ", m," : ",d['description'])
            i = i+1
        response = input("Choose one:")
        if response in choices:
            return response
        for m,d in choices.items():
            try:
                if d[None]==int(response):
                    return m
            except ValueError:
                pass
            
def chooseModelAndParameters(model):
    """
    Presents user with a choice of models and asks for parameters.
    """
    modelname = chooseModel(model)
    result = {None:modelname}
    model = Model(modelname)
    params = model.getParameterDictionary()
    del model
    for param in params.copy():
        if '.' in param:
            del params[param]
    for param,d in params.items():
        if d['type'][-4:]!="_Ptr":
            response = input(d['description']+" <"+ d['value'] + ">:")
            if response == '':
                result[param] = str(d['value'])
            else:
                result[param] = response
        else:
             result[param]= chooseModelandParameters(d['type'][:-4])
    return result

def createDictOfDict(theDict):
    """
    Creates a dictionary of dictionaries for a simple dictionary
    with dot-delimited keys. The result[gparent][parent][child]
    theDict["gparent.parent.child"]-->result["gparent"]["parent"]["child"]
    theDict["gparent.parent"]--> result[gparent][parent][None]
    """
    result = {}
    for p in theDict:
        if p==None:
            result[p] = theDict[p] 
        elif '.' in p:
            parent = p[:p.find('.')]
            child = p[p.find('.')+1:]
            if type(result[parent])!=dict:
                result[parent] = {None : result[parent]}
            result[parent][child] = theDict[p] 
        else:
            result[p] = theDict[p]
    for p in result:
        if type(result[p])==dict:
            result[p] = createDictOfDict(result[p])
    return result

def printDictOfDict(thedict):
    print(strDictOfDict(thedict))
    
def strDictOfDict(thedict,indent="",result=None):
    """
    Prints a dictionary of dictionaries in a readable format
    """
    if result==None:
        result = "{"
    else:
        result += "{"
    i=0
    for a in thedict:
        if i:
            result += ",\n"
        b = "'" + a + "'" if a!=None else "None"
        if type(thedict[a])==dict:
            if i:
                result += indent + b + " : "
            else:
                result += b + " : "
            result = strDictOfDict(thedict[a],indent+(len(a)+7)*' ',result)
        else:
            c = "'"+thedict[a]+"'"
            if i:
                result += indent + b + " : " + c
            else:
                result += b + " : " + c
        i = i+1
    
    result += "}"
    return result



def DialogGetModel(top):
    """
    Opens a dialog box with a selection of Models inheriting top and 
    returns the choice.

    Parameters
    ----------
    top : str
          The top level class from which to choose a child class

    Returns
    -------
    str : The name of the selected child class

    Requires
    --------
    tkinter
    """
    import tkinter
    from tkinter import ttk
    from tkinter import simpledialog

    class ModelDialog:
        # Opens a dialog box containing a bunch of choices. User 
        # double-clicks or presses enter to choose.
        # On completion, result contains the result.
        def ok(self,event):
            parameter = self.tr.selection()[0]
            self.result = parameter
            self.top.destroy()

        def __init__(self, title, choices):

            self.choices=choices

            top = self.top = tkinter.Tk()
            top.title(title)

            self.tr = tkinter.ttk.Treeview(top)
            self.tr["columns"] = ("Description")
            for t in self.tr["columns"]:
                self.tr.column(t)
                self.tr.heading(t,text=t)
                rows = 0
                cols = 0
            for c in self.choices:
                self.tr.insert("", "end", c, text=c,
                               values = (self.choices[c]["description"],))
                cols = max(cols,len(c)+len(self.choices[c]["description"]))
                rows += 1
                self.tr.pack(expand=True, fill="both")
                self.tr.bind('<Double-1>', self.ok)
                self.tr.bind('<Return>', self.ok)

            self.result = None
            top.geometry("%dx%d" % (cols*8,rows*22))
            top.wait_window(self.top)

    dialog = ModelDialog("Select model:", getModelDictionary(top))
    try:
        model = Model(dialog.result)
    except:
        model = DialogGetModel(top)

    return model
        

def DialogGetModelParameters(model):
    """
    Opens a dialog box with a selection of Models inheriting top and 
    returns the choice.

    Parameters
    ----------
    model : Model
            The model to fill parameters out with.
    
    Returns
    -------
    Model
    The model with new parameters

    Requires
    --------
    tkinter
    """
    import tkinter
    from tkinter import ttk
    from tkinter import simpledialog

    class ModelParameterDialog():
        # Creates a Treeview with all of a model's parameters. User can change 
        # values as desired. Upon completion, the model is set as desired.
        def __init__(self, model):
            self.model = model
            top = self.top = tkinter.Tk()
            self.top.title("Model Parameters")
            self.tree = tkinter.ttk.Treeview(self.top)
            self.tree["columns"] = ("value","description","type")
            self.tree.column("value")
            self.tree.column("description")
            self.tree.column("type")
            self.tree.heading("value", text="Value")
            self.tree.heading("description", text="Description")
            self.tree.heading("type", text="Type")
            self.reload()
            self.tree.pack(expand=True,fill='both')
            top.wait_window(self.top)

        def reload(self):
            x = self.tree.get_children()
            for t in x:
                self.tree.delete(t)
            params = self.model.getParameterDictionary()
            for p in params:
                if "_Ptr" in params[p]['type']:
                    if '.' in p:
                        parent = p[:p.rfind('.')]
                        child = p[p.rfind('.')+1:]
                    else:
                        parent = ""
                        child = p
                elif '.' in p:
                    parent = p[:p.rfind('.')]
                    child = p[p.rfind('.')+1:]
                else:
                    parent = ""
                    child = p

                self.tree.insert(parent,"end",p,text=child,
                                 values=(str(params[p]["value"]),
                                         params[p]["description"],
                                         params[p]["type"]))

            self.tree.bind('<Double-1>',self.clickon)
            self.tree.bind('<Return>',self.clickon)

        def clickon(self,event):
            parameter = self.tree.selection()[0]
            v = self.tree.item(parameter)['values']
            prompt = ("Parameter = " + parameter +
                      "\n" + str(v[1]) +
                      "\nDatatype = " + str(v[2]) +
                      "\nCurrent value = " + str(v[0]) +
                      "\nEnter new value:")
            if v[2][-4:]=='_Ptr':
                dialogresult = DialogGetModel(v[2][:-4])
                if (dialogresult!=None):
                    self.model.setParameters({parameter : dialogresult})
                    self.reload()
            else:
                answer = tkinter.simpledialog.askstring("Input", prompt,
                                                        parent=self.tree)
                if answer!="" and answer!=None:
                    temp = self.tree.item(parameter)['values']
                    temp[0]=answer
                    self.tree.item(parameter,values=temp)
                    self.model.setParameter(parameter,answer)
                    
    return ModelParameterDialog(model).model
  
