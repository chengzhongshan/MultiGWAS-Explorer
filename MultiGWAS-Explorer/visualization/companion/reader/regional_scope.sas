/* Caller context required by the unchanged original regional macros.
   Matches the retained reference submission's global declarations.
   CHECK_VAR_TYPE writes VAR_TYPE with PROC SQL INTO; a name first
   created inside that helper otherwise disappears on return. */
%global chr_var dotsize var_type;
