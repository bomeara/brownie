#-----------------------
# datatypes class(es)
# @ author Conrad Stack
#
# 'datatype' class will basically S4-wrap a character string taking 
# one of 4 or so values.
#
# datatype will be extended by two other classes: 'datatype_bound' and 'datatype_prob'
# These will eventually be used to represent confidence interval bounds ('datatype_bound') or 
# probabilities of certain states ('datatype_prob')
# 
#-----------------------


setClass("datatype",
		representation(type="character",name="character"),
		prototype=prototype(type=genericData(),name="NA"))

is_valid_datatype <- function(object)
{
	if(object@type %in% brownie.datatypes())
		return(TRUE)
	return(FALSE)
} 
setValidity("datatype",is_valid_datatype)


# Subtypes:
#----------

# adds level slot, which is the upper or lower density 
# level of the data. (e.g. 5%, 95%).  Clearly, this subtype 
# if for use with continuous data
#
setClass("datatype_bound",
		representation(level="numeric"),
		contains="datatype",
		prototype=prototype(level=100)
		)

is_valid_datatype_bound <- function(object)
{
	retval = TRUE
	if(object@type != contData())
		retval = FALSE
	
	if(object@level > 100.0 || object@level < 0.0)
		retval = FALSE
	
	return(retval)
} 
setValidity("datatype_bound",is_valid_datatype_bound)


# Adds state, prob, and ln to datatype.
# @prob is the probability that the @name variable
# is @state.  @ln indicates whether or not the probability
# is represented as the natural logarithm or not.
#
setClass("datatype_prob",
		representation(state="character",prob="numeric",ln="logical"),
		contains="datatype",
		prototype=prototype(ln=F))

is_valid_datatype_prob <- function(object)
{
	retval = TRUE
	
	if(object@type != discData())
		return(FALSE)
	
	# if prob is represented as log(prob)
	if(object@ln){
		if(object@prob > 0.0)
			retval = F
	} else {
		if(object@prob > 1.0 || object < 0.0)
			retval = F
	}
	
	return (retval)
}
setValidity("datatype_prob",is_valid_datatype_bound)

