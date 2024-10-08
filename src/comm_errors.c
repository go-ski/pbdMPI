/* Modified from "R-devel/main/src/error.c" and
 *               "R-devel/main/include/Rinternals.h".
 *
 * WCC: .Internal(...) is not allowed for users to modified and may change.
 *      However, some functions are exported in "Rinternals.h" which can be
 *      reused at somewhere. .External(...) can do the similar works,
 *      but the caller should be passed from R, since findCall() is not
 *      available in "Rinternals.h" and requires more R negative structures.
 *      Replace errorcall() by Rf_errorcall(), and warningcall() by
 *      Rf_warningcall().
 */

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

static int immediateWarning = 0;


/* Replace Rf_isValidString(). */
Rboolean api_R_isValidString(SEXP x){
    return ((TYPEOF(x) == STRSXP) && (LENGTH(x) > 0) && (TYPEOF(STRING_ELT(x, 0)) != NILSXP));
}

/* Replace Rf_errorcall() and Rf_warningcall(). */
// NORET void api_R_errorcall(SEXP, const char *, ...) R_PRINTF_FORMAT(2, 3);
// void api_R_warningcall(SEXP, const char *, ...) R_PRINTF_FORMAT(2, 3);


/* Origin: SEXP attribute_hidden do_stop(). */
SEXP api_R_stop(SEXP args){
	SEXP call, c_call;

	args = CDR(args);		/* get caller name */
	call = CAR(args);

	args = CDR(args);
	if(asLogical(CAR(args))){	/* find context -> "... in: ..:" */
		c_call = call;
	} else{
		c_call = R_NilValue;
	}

	args = CDR(args);
	if(CAR(args) != R_NilValue){	/* message */
		SETCAR(args, coerceVector(CAR(args), STRSXP));
		if(!api_R_isValidString(CAR(args))){
			errorcall(c_call,
				" [invalid string in comm.stop(.)]\n");
		}
		errorcall(c_call, "%s",
			translateChar(STRING_ELT(CAR(args), 0)));
	} else{
		errorcall(c_call, "\n");
	}

	return c_call;
} /* End of api_R_stop(). */


/* Origin: SEXP attribute_hidden do_warning(). */
SEXP api_R_warning(SEXP args){
	SEXP call, c_call;

	args = CDR(args);		/* get caller name */
	call = CAR(args);

	args = CDR(args);
	if(asLogical(CAR(args))){	/* find context -> "... in: ..:" */
		c_call = call;
	} else{
		c_call = R_NilValue;
	}

	args = CDR(args);

	if(asLogical(CAR(args))){	/* immediate = TRUE */
		immediateWarning = 1;
	} else{
		immediateWarning = 0;
	}

	args = CDR(args);
	if(CAR(args) != R_NilValue){
		SETCAR(args, coerceVector(CAR(args), STRSXP));
		if(!api_R_isValidString(CAR(args))){
			warningcall(c_call,
				" [invalid string in comm.warning(.)]\n");
		} else{
			warningcall(c_call, "%s",
				translateChar(STRING_ELT(CAR(args), 0)));
		}
	} else{
		warningcall(c_call, "%s", "");
	}
	immediateWarning = 0;	/* reset to internal calls */

	return CAR(args);
} /* End of api_R_warning(). */

