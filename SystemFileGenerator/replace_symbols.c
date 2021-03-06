/*
 * mex function which replaces symbols in a string. This function is used by the 
 * SystemFileGenerator. This mex function replaces replace_symbols_old.m. This 
 * mex function was written to speed up the generation of large systems.
 *
 * use: replace_symbols(<expression>, ...
 *                      <cell array containing char arrays of old symbols>, ...
 *                      <cell array containing char arrays of new symbols>)
 *
 * example:
 * replace_symbols( ...
 *      'sin(a * x * x * y)', {'x';'y';'z';'a'}, {'v_x', 'v_y', 'v_z', 'par_a'})
 *
 * returns 'sin(par_a * v_x * v_x * v_y)' 
 *
 * detailed specification:
 *
 * - A symbol that is to be replaced starts with a letter a-z or A-Z
 *
 * - A symbol contains only letters a-z, A-Z, underscores or digits 0-9.
 *
 * - A symbol in a string cannot be preceded by a letter, a digit, an 
 *   underscore, or a period. If a symbol is preceded by such a character, it is 
 *   ignored.
 *   The reason for this requirement is to prevent the replacement of 's' in 
 *   'cos(x)', if a user asks to replace 's'.
 *
 * - A symbol in a string cannot be followed by a letter, a digit or an 
 *   underscore. If a symbol is followed by such a character, then the symbol is 
 *   ignored.
 *   The reason for this requirement is to prevent the replacement of 's' in 
 *   'sin(x)', if a user asks to replace 's'.
 * 
 * - A symbol is only replaced once. Specifically replace_symbols scans the 
 *   string for all old symbols, and stores their locations. After scanning for
 *   old symbols, the detected old symbols are replaced by new symbols.
 *   for example: replace_symbols('a', {'a','b'}, {'b', 'c'}) returns 'b'
 * 
 */

#include <string.h>
#include <mex.h>
#include <stdbool.h>
#include <ctype.h>

#define STRING        0
#define OLD_SYMBOLS   1
#define NEW_SYMBOLS   2

typedef struct Node_tag {
  bool             is_symbol;
  int              si;
  char*            str;
  struct Node_tag* next;
} Node;


void check_cell(const mxArray* array, char* arrayname);
void check_1D  (const mxArray* array, char* arrayname);
void check_char(const mxArray* array, char* arrayname);
void check_null(const mxArray* array, char* arrayname);

void split_at_symbol(int s_i, Node* list);

char** old_symbols;
char** new_symbols;
int*   old_symbols_len;
int*   new_symbols_len;

void mexFunction(int n_output,       mxArray *mex_output[], 
                 int n_input,  const mxArray *mex_input []  ) {
  
  if (n_input != 3) {
    mexErrMsgIdAndTxt("replace_symbols:n_input",
                      "replace_symbols requires 3 arguments, but "
                      "%d are given. ", n_input);
  }
  
  if ( ! mxIsChar(mex_input[STRING]) ) {
    mexErrMsgIdAndTxt("replace_symbols:not_char", 
                      "The first argument is not a char array.");
  }
  
  if ( mxGetM(mex_input[STRING]) != 1 ) {
    mexErrMsgIdAndTxt("replace_symbols:first_arg_more_than_one_row", 
            "The first argument, "
            "which should be the string in which symbols are to be replaced, "
            "must have one row.");
  }
  
  char* input_string = (char*) mxArrayToString(mex_input[STRING]);
  
  char * old_symbols_description = "The second argument to replace_symbols, "
          "which should be a cell array containing the old symbols,";

  char * new_symbols_description = "The third argument to replace_symbols, "
          "which should be a cell array containing the new symbols,";
  
  check_cell(mex_input[OLD_SYMBOLS], old_symbols_description);
  check_1D  (mex_input[OLD_SYMBOLS], old_symbols_description);
  check_cell(mex_input[NEW_SYMBOLS], new_symbols_description);
  check_1D  (mex_input[NEW_SYMBOLS], new_symbols_description);
  
  int old_m = mxGetM(mex_input[OLD_SYMBOLS]);
  int old_n = mxGetN(mex_input[OLD_SYMBOLS]);
  int new_m = mxGetM(mex_input[NEW_SYMBOLS]);
  int new_n = mxGetN(mex_input[NEW_SYMBOLS]);
  
  if ( old_m * old_n != new_m * new_n ) {
    mexErrMsgIdAndTxt("replace_symbols:size_mismatch", 
            "The second argument and third arguments to replace_symbols, "
            "are not of the same size.");
  }
  
  int n_symbols = old_m * old_n;
  
  old_symbols     = mxMalloc(n_symbols * sizeof(char*));
  old_symbols_len = mxMalloc(n_symbols * sizeof(int  ));
  new_symbols     = mxMalloc(n_symbols * sizeof(char*));
  new_symbols_len = mxMalloc(n_symbols * sizeof(int  ));
  
  
  for (int i = 0; i < n_symbols; i++) {
    mxArray* old_symbol = mxGetCell(mex_input[OLD_SYMBOLS], i);
    mxArray* new_symbol = mxGetCell(mex_input[NEW_SYMBOLS], i);
    check_char(old_symbol, "One of the elements of the second input");
    check_char(new_symbol, "One of the elements of the third input");
    old_symbols    [i] = mxArrayToString(old_symbol);
    old_symbols_len[i] = mxGetN         (old_symbol);
    new_symbols    [i] = mxArrayToString(new_symbol);
    new_symbols_len[i] = mxGetN         (new_symbol);
  }
  
  Node* list = mxMalloc(sizeof(Node));
  list -> is_symbol = false;
  list -> str       = input_string;
  list -> next      = NULL;
  
  for (int i = 0; i < n_symbols; i++) {
    split_at_symbol(i, list);
  }
  
  int output_length = 0;
  
  {
    Node* node = list;
    do {
      if (node -> is_symbol) {
        output_length += new_symbols_len[node -> si];
      } else {
        output_length += strlen(node -> str);
      }
    } while (node = node -> next);
  }
  
  char* output_string = mxMalloc( (output_length + 1) * sizeof(char) );
  
  {
    Node* node = list;
    char* output_ptr = output_string;
    while(node) { 
      if (node -> is_symbol) {
        int si = node -> si;
        memcpy(output_ptr, new_symbols[si], new_symbols_len[si]);
        output_ptr += new_symbols_len[si];
      } else {
        memcpy(output_ptr, node -> str, strlen(node -> str));
        output_ptr += strlen(node -> str);
      }
      Node* next_node = node -> next;
      mxFree(node);
      node = next_node;
    }
    *output_ptr = '\0';
  }
  
  mex_output[0] = mxCreateString(output_string);

  for (int i = 0; i < n_symbols; i++ ) {
    mxFree(old_symbols[i]);
    mxFree(new_symbols[i]);
  }

  mxFree(old_symbols    );
  mxFree(old_symbols_len);
  mxFree(new_symbols    );
  mxFree(new_symbols_len);
  mxFree(input_string);
  mxFree(output_string);
  
}
  
void split_at_symbol(int si, Node* list) {
  Node* node = list;
  do {
    if (node -> is_symbol) {
      continue;
    }
    for (
            char* pch = strstr(node -> str, old_symbols[si]);
            pch;
            pch = strstr(pch + 1, old_symbols[si])) {
      bool symbol_continues_on_left;
      if (pch == list -> str) {
        symbol_continues_on_left = false;
      } else {
        symbol_continues_on_left = isalnum(pch[-1]) || pch[-1] == '_' 
                                                    || pch[-1] == '.';
        // pch[-1] == '.' is to prevent the 'e' possibly being replaced 
        // in 1.e-2 
      }
      if (symbol_continues_on_left) {
        continue;
      }
      
      bool symbol_continues_on_right;
      if (pch[old_symbols_len[si]] == '\0') {
        symbol_continues_on_right = false;
      } else {
        char c = pch[ old_symbols_len[si] ];
        symbol_continues_on_right = isalnum(c) || c == '_';
      }
      if (symbol_continues_on_right) {
        continue;
      }
      
    //  delta_len += new_symbols_len[si] - old_symbols_len[si];
      
      Node* next_in_list = node -> next;
      
      Node* suffix = mxMalloc(sizeof(Node));
      suffix -> is_symbol = false;
      suffix -> str       = pch + old_symbols_len[si];
      suffix -> next      = next_in_list;
      
      Node* symbol = mxMalloc(sizeof(Node));
      symbol -> is_symbol = true;
      symbol -> si        = si;
      symbol -> next      = suffix;
      
      node   -> next = symbol;
      node = suffix;
      *pch++ = '\0';
    } // end of for ( pch ...
  } while ( node = node -> next );
}

void check_cell(const mxArray* array, char* arrayname) {
  check_null(array, arrayname);
  if ( ! mxIsCell(array) ) {
    mexErrMsgIdAndTxt("replace_symbols:not_cell", 
                      "%s is not a cell array.", arrayname);
  }
}
    
void check_1D(const mxArray* array, char* arrayname) {
  check_null(array, arrayname);
  if ( ! (mxGetM(array) == 1 || mxGetN(array) == 1) ) {
    mexErrMsgIdAndTxt("replace_symbols:not_1D", 
                      "%s is not a one-dimensional array.", arrayname);
  }
}

void check_char(const mxArray* array, char* arrayname) {
  check_null(array, arrayname);
  if ( ! mxIsChar(array) ) {
    mexErrMsgIdAndTxt("replace_symbols:not_char", 
                      "%s is not a char array.", arrayname);
  }
}

void check_null(const mxArray* array, char* arrayname) {
  if ( ! array ) {
    mexErrMsgIdAndTxt("replace_symbols:null", 
                      "%s is NULL (i.e. empty)", arrayname);
  }
}
