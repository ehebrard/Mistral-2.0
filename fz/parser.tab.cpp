/* A Bison parser, made by GNU Bison 3.0.4.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "3.0.4"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 1

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* Copy the first part of user declarations.  */
#line 40 "parser.yxx" /* yacc.c:339  */

#define YYPARSE_PARAM parm
#define YYLEX_PARAM static_cast<ParserState*>(parm)->yyscanner
#include "flatzinc.hpp"
#include "parser.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>

#ifdef HAVE_MMAP
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#endif

using namespace std;

int yyparse(void*);
int yylex(YYSTYPE*, void* scanner);
int yylex_init (void** scanner);
int yylex_destroy (void* scanner);
int yyget_lineno (void* scanner);
void yyset_extra (void* user_defined ,void* yyscanner );

extern int yydebug;

using namespace FlatZinc;

void yyerror(void* parm, const char *str) {
  ParserState* pp = static_cast<ParserState*>(parm);
  pp->err << "Error: " << str
          << " in line no. " << yyget_lineno(pp->yyscanner)
          << std::endl;
  pp->hadError = true;
}

void yyassert(ParserState* pp, bool cond, const char* str)
{
  if (!cond) {
    pp->err << "Error: " << str
            << " in line no. " << yyget_lineno(pp->yyscanner)
            << std::endl;
    pp->hadError = true;
  }
}

/*
 * The symbol tables
 *
 */

AST::Node* getArrayElement(ParserState* pp, string id, unsigned int offset) {
  if (offset > 0) {
    vector<int> tmp;
    if (pp->intvararrays.get(id, tmp) && offset<=tmp.size())
      return new AST::IntVar(tmp[offset-1]);
    if (pp->boolvararrays.get(id, tmp) && offset<=tmp.size())
      return new AST::BoolVar(tmp[offset-1]);
    if (pp->setvararrays.get(id, tmp) && offset<=tmp.size())
      return new AST::SetVar(tmp[offset-1]);

    if (pp->intvalarrays.get(id, tmp) && offset<=tmp.size())
      return new AST::IntLit(tmp[offset-1]);
    if (pp->boolvalarrays.get(id, tmp) && offset<=tmp.size())
      return new AST::BoolLit(tmp[offset-1]);
    vector<AST::SetLit> tmpS;
    if (pp->setvalarrays.get(id, tmpS) && offset<=tmpS.size())
      return new AST::SetLit(tmpS[offset-1]);
  }

  pp->err << "Error: array access to " << id << " invalid"
          << " in line no. "
          << yyget_lineno(pp->yyscanner) << std::endl;
  pp->hadError = true;
  return new AST::IntVar(0); // keep things consistent
}
AST::Node* getVarRefArg(ParserState* pp, string id, bool annotation = false) {
  int tmp;
  if (pp->intvarTable.get(id, tmp))
    return new AST::IntVar(tmp);
  if (pp->boolvarTable.get(id, tmp))
    return new AST::BoolVar(tmp);
  if (pp->setvarTable.get(id, tmp))
    return new AST::SetVar(tmp);
  if (annotation)
    return new AST::Atom(id);
  pp->err << "Error: undefined variable " << id
          << " in line no. "
          << yyget_lineno(pp->yyscanner) << std::endl;
  pp->hadError = true;
  return new AST::IntVar(0); // keep things consistent
}

void addDomainConstraint(ParserState* pp, std::string id, AST::Node* var,
                         Option<AST::SetLit* >& dom) {
  if (!dom())
    return;
  AST::Array* args = new AST::Array(2);
  args->a[0] = var;
  args->a[1] = dom.some();
  pp->domainConstraints.push_back(new ConExpr(id, args));
}

/*
 * Initialize the root gecode space
 *
 */

void initfg(ParserState* pp) {
  if (!pp->hadError)
    pp->fg->init(pp->intvars.size(),
                 pp->boolvars.size(),
                 pp->setvars.size());

  for (unsigned int i=0; i<pp->intvars.size(); i++) {
    if (!pp->hadError) {
      try {
        pp->fg->newIntVar(static_cast<IntVarSpec*>(pp->intvars[i].second));
      } catch (FlatZinc::Error& e) {
        yyerror(pp, e.toString().c_str());
      }
    }
    if (pp->intvars[i].first[0] != '[') {
      delete pp->intvars[i].second;
      pp->intvars[i].second = NULL;
    }
  }
  for (unsigned int i=0; i<pp->boolvars.size(); i++) {
    if (!pp->hadError) {
      try {
        pp->fg->newBoolVar(
          static_cast<BoolVarSpec*>(pp->boolvars[i].second));
      } catch (FlatZinc::Error& e) {
        yyerror(pp, e.toString().c_str());
      }
    }
    if (pp->boolvars[i].first[0] != '[') {
      delete pp->boolvars[i].second;
      pp->boolvars[i].second = NULL;
    }
  }
  for (unsigned int i=0; i<pp->setvars.size(); i++) {
    if (!pp->hadError) {
      try {
        pp->fg->newSetVar(static_cast<SetVarSpec*>(pp->setvars[i].second));
      } catch (FlatZinc::Error& e) {
        yyerror(pp, e.toString().c_str());
      }
    }
    if (pp->setvars[i].first[0] != '[') {
      delete pp->setvars[i].second;
      pp->setvars[i].second = NULL;
    }
  }
  for (unsigned int i=pp->domainConstraints.size(); i--;) {
    if (!pp->hadError) {
      try {
        assert(pp->domainConstraints[i]->args->a.size() == 2);
        pp->fg->postConstraint(*pp->domainConstraints[i], NULL);
        delete pp->domainConstraints[i];
      } catch (FlatZinc::Error& e) {
        yyerror(pp, e.toString().c_str());
      }
    }
  }
}

void fillPrinter(ParserState& pp, FlatZinc::Printer& p) {
  p.init(pp.getOutput());
}

AST::Node* arrayOutput(AST::Call* ann) {
  AST::Array* a = NULL;

  if (ann->args->isArray()) {
    a = ann->args->getArray();
  } else {
    a = new AST::Array(ann->args);
  }

  std::ostringstream oss;

  oss << "array" << a->a.size() << "d(";
  for (unsigned int i=0; i<a->a.size(); i++) {
    AST::SetLit* s = a->a[i]->getSet();
    if (s->empty())
      oss << "{}, ";
    else if (s->interval)
      oss << s->min << ".." << s->max << ", ";
    else {
      oss << "{";
      for (unsigned int j=0; j<s->s.size(); j++) {
        oss << s->s[j];
        if (j<s->s.size()-1)
          oss << ",";
      }
      oss << "}, ";
    }
  }

  if (!ann->args->isArray()) {
    a->a[0] = NULL;
    delete a;
  }
  return new AST::String(oss.str());
}

/*
 * The main program
 *
 */

namespace FlatZinc {

  FlatZincModel* parse(const std::string& filename,
                       Solver& solver,
                       Printer& p, std::ostream& err,
                       FlatZincModel* fzs) {
#ifdef HAVE_MMAP
    int fd;
    char* data;
    struct stat sbuf;
    fd = open(filename.c_str(), O_RDONLY);
    if (fd == -1) {
      err << "Cannot open file " << filename << endl;
      return NULL;
    }
    if (stat(filename.c_str(), &sbuf) == -1) {
      err << "Cannot stat file " << filename << endl;
      return NULL;
    }
    data = (char*)mmap((caddr_t)0, sbuf.st_size, PROT_READ, MAP_SHARED, fd,0);
    if (data == (caddr_t)(-1)) {
      err << "Cannot mmap file " << filename << endl;
      return NULL;
    }

    if (fzs == NULL) {
      fzs = new FlatZincModel();
    }
    ParserState pp(data, sbuf.st_size, err, fzs);
#else
    std::ifstream file;
    file.open(filename.c_str());
    if (!file.is_open()) {
      err << "Cannot open file " << filename << endl;
      return NULL;
    }
    std::string s = string(istreambuf_iterator<char>(file),
                           istreambuf_iterator<char>());
    if (fzs == NULL) {
      fzs = new FlatZincModel(solver);
    }
    ParserState pp(s, err, fzs);
#endif
    yylex_init(&pp.yyscanner);
    yyset_extra(&pp, pp.yyscanner);
    // yydebug = 1;
    yyparse(&pp);
    fillPrinter(pp, p);

    if (pp.yyscanner)
      yylex_destroy(pp.yyscanner);
    return pp.hadError ? NULL : pp.fg;
  }

  FlatZincModel* parse(std::istream& is,
                       Solver& solver,
                       Printer& p, std::ostream& err,
                       FlatZincModel* fzs) {
    std::string s = string(istreambuf_iterator<char>(is),
                           istreambuf_iterator<char>());

    if (fzs == NULL) {
      fzs = new FlatZincModel(solver);
    }
    ParserState pp(s, err, fzs);
    yylex_init(&pp.yyscanner);
    yyset_extra(&pp, pp.yyscanner);
    // yydebug = 1;
    yyparse(&pp);
    fillPrinter(pp, p);

    if (pp.yyscanner)
      yylex_destroy(pp.yyscanner);
    return pp.hadError ? NULL : pp.fg;
  }

}


#line 363 "parser.tab.cpp" /* yacc.c:339  */

# ifndef YY_NULLPTR
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULLPTR nullptr
#  else
#   define YY_NULLPTR 0
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 1
#endif

/* In a future release of Bison, this section will be replaced
   by #include "parser.tab.hpp".  */
#ifndef YY_YY_PARSER_TAB_HPP_INCLUDED
# define YY_YY_PARSER_TAB_HPP_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 1
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    FZ_INT_LIT = 258,
    FZ_BOOL_LIT = 259,
    FZ_FLOAT_LIT = 260,
    FZ_ID = 261,
    FZ_STRING_LIT = 262,
    FZ_VAR = 263,
    FZ_PAR = 264,
    FZ_ANNOTATION = 265,
    FZ_ANY = 266,
    FZ_ARRAY = 267,
    FZ_BOOL = 268,
    FZ_CASE = 269,
    FZ_COLONCOLON = 270,
    FZ_CONSTRAINT = 271,
    FZ_DEFAULT = 272,
    FZ_DOTDOT = 273,
    FZ_ELSE = 274,
    FZ_ELSEIF = 275,
    FZ_ENDIF = 276,
    FZ_ENUM = 277,
    FZ_FLOAT = 278,
    FZ_FUNCTION = 279,
    FZ_IF = 280,
    FZ_INCLUDE = 281,
    FZ_INT = 282,
    FZ_LET = 283,
    FZ_MAXIMIZE = 284,
    FZ_MINIMIZE = 285,
    FZ_OF = 286,
    FZ_SATISFY = 287,
    FZ_OUTPUT = 288,
    FZ_PREDICATE = 289,
    FZ_RECORD = 290,
    FZ_SET = 291,
    FZ_SHOW = 292,
    FZ_SHOWCOND = 293,
    FZ_SOLVE = 294,
    FZ_STRING = 295,
    FZ_TEST = 296,
    FZ_THEN = 297,
    FZ_TUPLE = 298,
    FZ_TYPE = 299,
    FZ_VARIANT_RECORD = 300,
    FZ_WHERE = 301
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 337 "parser.yxx" /* yacc.c:355  */
 int iValue; char* sValue; bool bValue; double dValue;
         std::vector<int>* setValue;
         FlatZinc::AST::SetLit* setLit;
         std::vector<double>* floatSetValue;
         std::vector<FlatZinc::AST::SetLit>* setValueList;
         FlatZinc::Option<FlatZinc::AST::SetLit* > oSet;
         FlatZinc::VarSpec* varSpec;
         FlatZinc::Option<FlatZinc::AST::Node*> oArg;
         std::vector<FlatZinc::VarSpec*>* varSpecVec;
         FlatZinc::Option<std::vector<FlatZinc::VarSpec*>* > oVarSpecVec;
         FlatZinc::AST::Node* arg;
         FlatZinc::AST::Array* argVec;
       

#line 465 "parser.tab.cpp" /* yacc.c:355  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif



int yyparse (void *parm);

#endif /* !YY_YY_PARSER_TAB_HPP_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 481 "parser.tab.cpp" /* yacc.c:358  */

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

#ifndef YY_ATTRIBUTE
# if (defined __GNUC__                                               \
      && (2 < __GNUC__ || (__GNUC__ == 2 && 96 <= __GNUC_MINOR__)))  \
     || defined __SUNPRO_C && 0x5110 <= __SUNPRO_C
#  define YY_ATTRIBUTE(Spec) __attribute__(Spec)
# else
#  define YY_ATTRIBUTE(Spec) /* empty */
# endif
#endif

#ifndef YY_ATTRIBUTE_PURE
# define YY_ATTRIBUTE_PURE   YY_ATTRIBUTE ((__pure__))
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# define YY_ATTRIBUTE_UNUSED YY_ATTRIBUTE ((__unused__))
#endif

#if !defined _Noreturn \
     && (!defined __STDC_VERSION__ || __STDC_VERSION__ < 201112)
# if defined _MSC_VER && 1200 <= _MSC_VER
#  define _Noreturn __declspec (noreturn)
# else
#  define _Noreturn YY_ATTRIBUTE ((__noreturn__))
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

#if defined __GNUC__ && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN \
    _Pragma ("GCC diagnostic push") \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")\
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# define YY_IGNORE_MAYBE_UNINITIALIZED_END \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif


#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYSIZE_T yynewbytes;                                            \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / sizeof (*yyptr);                          \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  7
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   324

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  57
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  66
/* YYNRULES -- Number of rules.  */
#define YYNRULES  154
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  334

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   301

#define YYTRANSLATE(YYX)                                                \
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, without out-of-bounds checking.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      48,    49,     2,     2,    50,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    51,    47,
       2,    54,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    52,     2,    53,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    55,     2,    56,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   436,   436,   438,   440,   443,   444,   448,   449,   453,
     454,   456,   458,   461,   462,   469,   471,   473,   476,   477,
     480,   483,   484,   485,   486,   489,   490,   491,   492,   495,
     496,   499,   500,   507,   537,   566,   571,   601,   625,   634,
     646,   705,   757,   764,   819,   832,   845,   852,   866,   870,
     885,   909,   910,   914,   916,   919,   919,   921,   925,   927,
     942,   966,   967,   971,   973,   977,   981,   983,   998,  1022,
    1023,  1027,  1029,  1032,  1035,  1037,  1052,  1076,  1077,  1081,
    1083,  1086,  1091,  1092,  1097,  1098,  1103,  1104,  1109,  1110,
    1114,  1128,  1141,  1163,  1165,  1167,  1173,  1175,  1188,  1189,
    1196,  1198,  1205,  1206,  1210,  1212,  1217,  1218,  1222,  1224,
    1229,  1230,  1234,  1236,  1241,  1242,  1246,  1248,  1256,  1258,
    1262,  1264,  1269,  1270,  1274,  1276,  1278,  1280,  1282,  1331,
    1345,  1346,  1350,  1352,  1360,  1371,  1393,  1394,  1402,  1403,
    1407,  1409,  1413,  1417,  1421,  1423,  1427,  1429,  1433,  1435,
    1437,  1439,  1441,  1484,  1495
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 1
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "FZ_INT_LIT", "FZ_BOOL_LIT",
  "FZ_FLOAT_LIT", "FZ_ID", "FZ_STRING_LIT", "FZ_VAR", "FZ_PAR",
  "FZ_ANNOTATION", "FZ_ANY", "FZ_ARRAY", "FZ_BOOL", "FZ_CASE",
  "FZ_COLONCOLON", "FZ_CONSTRAINT", "FZ_DEFAULT", "FZ_DOTDOT", "FZ_ELSE",
  "FZ_ELSEIF", "FZ_ENDIF", "FZ_ENUM", "FZ_FLOAT", "FZ_FUNCTION", "FZ_IF",
  "FZ_INCLUDE", "FZ_INT", "FZ_LET", "FZ_MAXIMIZE", "FZ_MINIMIZE", "FZ_OF",
  "FZ_SATISFY", "FZ_OUTPUT", "FZ_PREDICATE", "FZ_RECORD", "FZ_SET",
  "FZ_SHOW", "FZ_SHOWCOND", "FZ_SOLVE", "FZ_STRING", "FZ_TEST", "FZ_THEN",
  "FZ_TUPLE", "FZ_TYPE", "FZ_VARIANT_RECORD", "FZ_WHERE", "';'", "'('",
  "')'", "','", "':'", "'['", "']'", "'='", "'{'", "'}'", "$accept",
  "model", "preddecl_items", "preddecl_items_head", "vardecl_items",
  "vardecl_items_head", "constraint_items", "constraint_items_head",
  "preddecl_item", "pred_arg_list", "pred_arg_list_head", "pred_arg",
  "pred_arg_type", "pred_arg_simple_type", "pred_array_init",
  "pred_array_init_arg", "vardecl_item", "int_init", "int_init_list",
  "int_init_list_head", "list_tail", "int_var_array_literal", "float_init",
  "float_init_list", "float_init_list_head", "float_var_array_literal",
  "bool_init", "bool_init_list", "bool_init_list_head",
  "bool_var_array_literal", "set_init", "set_init_list",
  "set_init_list_head", "set_var_array_literal",
  "vardecl_int_var_array_init", "vardecl_bool_var_array_init",
  "vardecl_float_var_array_init", "vardecl_set_var_array_init",
  "constraint_item", "solve_item", "int_ti_expr_tail", "bool_ti_expr_tail",
  "float_ti_expr_tail", "set_literal", "int_list", "int_list_head",
  "bool_list", "bool_list_head", "float_list", "float_list_head",
  "set_literal_list", "set_literal_list_head", "flat_expr_list",
  "flat_expr", "non_array_expr_opt", "non_array_expr",
  "non_array_expr_list", "non_array_expr_list_head", "solve_expr",
  "minmax", "annotations", "annotations_head", "annotation",
  "annotation_list", "annotation_expr", "ann_non_array_expr", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,    59,    40,    41,
      44,    58,    91,    93,    61,   123,   125
};
# endif

#define YYPACT_NINF -122

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-122)))

#define YYTABLE_NINF -1

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
     -17,    18,    48,    30,   -17,    28,    36,  -122,    63,    90,
      35,    49,  -122,    71,   109,    98,    30,    69,    81,    73,
    -122,    27,   131,  -122,  -122,   105,   146,    89,    92,    93,
     138,   140,    33,  -122,    97,   108,   141,   125,    98,   116,
     121,  -122,   159,  -122,   112,   117,  -122,  -122,   145,   122,
     127,  -122,   133,  -122,  -122,  -122,    33,  -122,  -122,   135,
     136,   181,   182,   183,   172,   176,   142,  -122,   189,  -122,
     147,   176,   149,   151,  -122,  -122,   176,  -122,    19,    33,
    -122,    27,  -122,   194,   152,   198,   148,   200,   150,   176,
     176,   176,   204,    22,   154,   195,   203,  -122,   118,    75,
    -122,  -122,   157,   196,  -122,    -6,  -122,  -122,  -122,  -122,
     206,  -122,  -122,  -122,  -122,   163,   163,   163,   165,   202,
    -122,  -122,   -34,  -122,    22,   109,  -122,  -122,  -122,  -122,
     177,    22,   176,   202,  -122,  -122,   169,   177,  -122,    21,
    -122,  -122,  -122,  -122,  -122,   216,   177,   220,    19,   197,
     176,   177,  -122,  -122,  -122,   199,   222,    22,   124,  -122,
       3,   170,  -122,  -122,   173,   177,  -122,   180,   186,   176,
     118,   185,  -122,  -122,  -122,  -122,    56,   163,  -122,    83,
    -122,    40,   187,   188,    22,  -122,  -122,   177,   190,  -122,
     177,  -122,  -122,  -122,   235,   112,  -122,  -122,   139,   191,
     193,   209,   201,  -122,  -122,  -122,  -122,  -122,  -122,   192,
    -122,   215,   205,   208,   210,   241,   242,    33,   243,  -122,
      33,   244,   245,   247,   176,   176,   212,   176,   213,   176,
     176,   176,   211,   214,   248,   217,   249,   218,   219,   221,
     224,   225,   176,   226,   176,   227,  -122,   228,  -122,   229,
    -122,   253,   255,   230,   109,   231,    55,  -122,    91,  -122,
     103,  -122,   233,   135,   234,   136,   236,   237,   239,  -122,
    -122,   240,  -122,   246,   232,  -122,   250,  -122,   251,   256,
    -122,   257,  -122,   252,   258,  -122,  -122,  -122,  -122,    10,
    -122,    46,  -122,   263,  -122,    55,  -122,   264,  -122,    91,
    -122,   266,  -122,   103,  -122,   202,  -122,   254,   260,   259,
    -122,   261,   262,  -122,   265,  -122,   267,  -122,   268,  -122,
    -122,    10,  -122,   271,  -122,    46,  -122,  -122,  -122,  -122,
    -122,   269,  -122,  -122
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       3,     0,     0,     7,     4,     0,     0,     1,     0,     0,
       0,     0,    93,     0,   102,    11,     8,     0,     0,     0,
       5,    16,     0,    96,    98,     0,   102,     0,     0,     0,
       0,     0,     0,   104,     0,    55,     0,     0,    12,     0,
       0,     9,     0,     6,     0,     0,    27,    28,     0,     0,
      55,    18,     0,    24,    25,    95,     0,   108,   112,    55,
      55,     0,     0,     0,     0,   138,     0,    94,    56,   103,
       0,   138,     0,     0,    13,    10,   138,    23,     0,     0,
      15,    56,    17,     0,     0,    56,     0,    56,     0,   138,
     138,   138,     0,     0,     0,   139,     0,   105,     0,     0,
       2,    14,     0,     0,    31,     0,    29,    26,    19,    20,
       0,   109,    97,   113,    99,   122,   122,   122,     0,   149,
     148,   150,   152,   154,     0,   102,   151,   140,   143,   146,
       0,     0,   138,   125,   124,   126,   128,   130,   127,     0,
     118,   120,   137,   136,    91,     0,     0,     0,     0,     0,
     138,     0,    33,    34,    35,     0,     0,     0,     0,   144,
       0,     0,    38,   141,     0,     0,   132,     0,    55,   138,
       0,   134,    92,    37,    32,    30,     0,   122,   123,     0,
     101,     0,   152,     0,     0,   147,   100,     0,     0,   121,
      56,   131,    90,   119,     0,     0,    21,    36,     0,     0,
       0,     0,     0,   142,   153,   145,    39,   129,   133,     0,
      22,     0,     0,     0,     0,     0,     0,     0,     0,   135,
       0,     0,     0,     0,   138,   138,     0,   138,     0,   138,
     138,   138,     0,     0,     0,     0,     0,    82,    84,    86,
       0,     0,   138,     0,   138,     0,    40,     0,    41,     0,
      42,   106,   110,     0,   102,    88,    51,    83,    69,    85,
      61,    87,     0,    55,     0,    55,     0,     0,     0,    43,
      48,    49,    53,     0,    55,    66,    67,    71,     0,    55,
      58,    59,    63,     0,    55,    45,   107,    46,   111,   114,
      44,    77,    89,     0,    57,    56,    52,     0,    73,    56,
      70,     0,    65,    56,    62,     0,   116,     0,    55,    75,
      79,     0,    55,    74,     0,    54,     0,    72,     0,    64,
      47,    56,   115,     0,    81,    56,    78,    50,    68,    60,
     117,     0,    80,    76
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -122,  -122,  -122,  -122,  -122,  -122,  -122,  -122,   279,  -122,
    -122,   238,  -122,   -39,  -122,   153,   273,   -25,  -122,  -122,
     -50,  -122,   -10,  -122,  -122,  -122,    -5,  -122,  -122,  -122,
     -30,  -122,  -122,  -122,  -122,  -122,  -122,  -122,   275,  -122,
      -1,    99,   100,   -90,  -121,  -122,  -122,    45,  -122,    51,
    -122,  -122,  -122,   130,  -105,  -114,  -122,  -122,  -122,  -122,
     -70,  -122,   -86,   158,  -122,   166
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,     2,     3,     4,    15,    16,    37,    38,     5,    49,
      50,    51,    52,    53,   105,   106,    17,   272,   273,   274,
      69,   257,   282,   283,   284,   261,   277,   278,   279,   259,
     310,   311,   312,   292,   246,   248,   250,   269,    39,    72,
      54,    28,    29,   138,    34,    35,   262,    59,   264,    60,
     307,   308,   139,   140,   152,   141,   167,   168,   172,   145,
      94,    95,   159,   160,   128,   129
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_uint16 yytable[] =
{
      82,    99,    18,   126,   161,    77,   102,   127,    27,    86,
      88,   153,   154,   305,   157,    18,   162,     1,   158,   115,
     116,   117,   103,   166,     6,   119,   120,   121,   122,   123,
       8,    66,   173,     8,   126,    44,     8,   178,     9,    45,
      46,   126,    10,    11,   148,   163,   104,   149,     7,   305,
      47,   188,   309,   184,    12,    84,   185,    12,   270,     8,
      12,   271,   164,    48,   195,   125,    13,   126,   126,    46,
     169,   170,   197,   206,   124,    20,   208,   125,   107,    47,
     177,    22,    14,    12,    21,    14,     8,    30,    14,   203,
     184,   198,    48,     8,   126,   275,   199,   276,   205,   192,
      31,   125,    32,    23,   142,   143,   200,   144,   280,   281,
      12,    14,    33,    24,    36,     8,    41,    12,   191,   201,
      43,   133,   134,   135,   136,    46,    25,   119,   120,   121,
     182,   123,    42,   267,    55,    47,    56,   196,    14,    12,
      61,    64,     8,    62,    63,    26,    65,    70,    48,    33,
      57,    58,    23,    67,   232,   233,   210,   235,    68,   237,
     238,   239,    24,    74,    71,    76,    12,    14,    75,    78,
     137,    80,   253,   125,   255,   211,    79,    81,   202,   125,
     133,   134,   135,   136,    83,    85,    87,    89,    90,    91,
      92,    93,    97,    96,    26,    98,   100,   212,   101,   306,
     109,   313,   111,   110,   112,   113,   114,   118,   130,   132,
     131,   146,   150,   286,   147,   288,   226,   151,   155,   228,
     156,   165,   171,   174,   296,   180,   186,   187,   176,   300,
     179,   330,   125,   189,   304,   313,   190,   194,   209,   158,
     217,   204,   215,   207,   216,   219,   220,   224,   225,   227,
     229,   230,   218,   231,   242,   244,   221,    57,   322,   222,
      58,   223,   326,   234,   236,   240,   314,   316,   241,   318,
     315,   243,   245,   247,   331,   249,   251,   252,   254,   256,
     258,   260,   295,    19,   266,   268,   285,   287,   289,    40,
     290,   291,   293,   319,   317,   332,   263,   213,   214,   294,
     193,   175,   297,   265,   298,   302,   299,   320,   303,   301,
     321,   323,   325,    73,   324,   181,     0,     0,   327,   108,
     328,   329,   333,     0,   183
};

static const yytype_int16 yycheck[] =
{
      50,    71,     3,    93,   125,    44,    76,    93,     9,    59,
      60,   116,   117,     3,    48,    16,   130,    34,    52,    89,
      90,    91,     3,   137,     6,     3,     4,     5,     6,     7,
       3,    32,   146,     3,   124,     8,     3,   151,     8,    12,
      13,   131,    12,    13,    50,   131,    27,    53,     0,     3,
      23,   165,     6,    50,    27,    56,    53,    27,     3,     3,
      27,     6,   132,    36,     8,    55,    36,   157,   158,    13,
      49,    50,   177,   187,    52,    47,   190,    55,    79,    23,
     150,    18,    55,    27,    48,    55,     3,    52,    55,    49,
      50,     8,    36,     3,   184,     4,    13,     6,   184,   169,
      51,    55,    31,    13,    29,    30,    23,    32,     5,     6,
      27,    55,     3,    23,    16,     3,    47,    27,   168,    36,
      47,     3,     4,     5,     6,    13,    36,     3,     4,     5,
       6,     7,    51,   254,     3,    23,    31,   176,    55,    27,
      51,     3,     3,    51,    51,    55,     6,     6,    36,     3,
       4,     5,    13,    56,   224,   225,   195,   227,    50,   229,
     230,   231,    23,    47,    39,     6,    27,    55,    47,    52,
      52,    49,   242,    55,   244,    36,    31,    50,   179,    55,
       3,     4,     5,     6,    51,    50,    50,     6,     6,     6,
      18,    15,     3,    51,    55,    48,    47,   198,    47,   289,
       6,   291,     4,    51,    56,     5,    56,     3,    54,     6,
      15,    54,     6,   263,    18,   265,   217,    54,    53,   220,
      18,    52,     6,     3,   274,     3,    56,    54,    31,   279,
      31,   321,    55,    53,   284,   325,    50,    52,     3,    52,
      31,    53,    51,    53,    51,    53,    31,     6,     6,     6,
       6,     6,    51,     6,     6,     6,    51,     4,   308,    51,
       5,    51,   312,    51,    51,    54,     3,     3,    54,     3,
     295,    54,    54,    54,     3,    54,    52,    52,    52,    52,
      52,    52,    50,     4,    54,    54,    53,    53,    52,    16,
      53,    52,    52,   303,   299,   325,   251,   198,   198,    53,
     170,   148,    52,   252,    53,    53,    50,    53,    50,    52,
      50,    52,    50,    38,    53,   157,    -1,    -1,    53,    81,
      53,    53,    53,    -1,   158
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,    34,    58,    59,    60,    65,     6,     0,     3,     8,
      12,    13,    27,    36,    55,    61,    62,    73,    97,    65,
      47,    48,    18,    13,    23,    36,    55,    97,    98,    99,
      52,    51,    31,     3,   101,   102,    16,    63,    64,    95,
      73,    47,    51,    47,     8,    12,    13,    23,    36,    66,
      67,    68,    69,    70,    97,     3,    31,     4,     5,   104,
     106,    51,    51,    51,     3,     6,    97,    56,    50,    77,
       6,    39,    96,    95,    47,    47,     6,    70,    52,    31,
      49,    50,    77,    51,    97,    50,    77,    50,    77,     6,
       6,     6,    18,    15,   117,   118,    51,     3,    48,   117,
      47,    47,   117,     3,    27,    71,    72,    97,    68,     6,
      51,     4,    56,     5,    56,   117,   117,   117,     3,     3,
       4,     5,     6,     7,    52,    55,   100,   119,   121,   122,
      54,    15,     6,     3,     4,     5,     6,    52,   100,   109,
     110,   112,    29,    30,    32,   116,    54,    18,    50,    53,
       6,    54,   111,   111,   111,    53,    18,    48,    52,   119,
     120,   101,   112,   119,   117,    52,   112,   113,   114,    49,
      50,     6,   115,   112,     3,    72,    31,   117,   112,    31,
       3,   120,     6,   122,    50,    53,    56,    54,   112,    53,
      50,    77,   117,   110,    52,     8,    70,   111,     8,    13,
      23,    36,    97,    49,    53,   119,   112,    53,   112,     3,
      70,    36,    97,    98,    99,    51,    51,    31,    51,    53,
      31,    51,    51,    51,     6,     6,    97,     6,    97,     6,
       6,     6,   117,   117,    51,   117,    51,   117,   117,   117,
      54,    54,     6,    54,     6,    54,    91,    54,    92,    54,
      93,    52,    52,   117,    52,   117,    52,    78,    52,    86,
      52,    82,   103,   104,   105,   106,    54,   101,    54,    94,
       3,     6,    74,    75,    76,     4,     6,    83,    84,    85,
       5,     6,    79,    80,    81,    53,    77,    53,    77,    52,
      53,    52,    90,    52,    53,    50,    77,    52,    53,    50,
      77,    52,    53,    50,    77,     3,   100,   107,   108,     6,
      87,    88,    89,   100,     3,    74,     3,    83,     3,    79,
      53,    50,    77,    52,    53,    50,    77,    53,    53,    53,
     100,     3,    87,    53
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    57,    58,    59,    59,    60,    60,    61,    61,    62,
      62,    63,    63,    64,    64,    65,    66,    66,    67,    67,
      68,    69,    69,    69,    69,    70,    70,    70,    70,    71,
      71,    72,    72,    73,    73,    73,    73,    73,    73,    73,
      73,    73,    73,    73,    73,    73,    73,    73,    74,    74,
      74,    75,    75,    76,    76,    77,    77,    78,    79,    79,
      79,    80,    80,    81,    81,    82,    83,    83,    83,    84,
      84,    85,    85,    86,    87,    87,    87,    88,    88,    89,
      89,    90,    91,    91,    92,    92,    93,    93,    94,    94,
      95,    96,    96,    97,    97,    97,    98,    98,    99,    99,
     100,   100,   101,   101,   102,   102,   103,   103,   104,   104,
     105,   105,   106,   106,   107,   107,   108,   108,   109,   109,
     110,   110,   111,   111,   112,   112,   112,   112,   112,   112,
     113,   113,   114,   114,   115,   115,   116,   116,   117,   117,
     118,   118,   119,   119,   120,   120,   121,   121,   122,   122,
     122,   122,   122,   122,   122
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     5,     0,     1,     2,     3,     0,     1,     2,
       3,     0,     1,     2,     3,     5,     0,     2,     1,     3,
       3,     6,     7,     2,     1,     1,     3,     1,     1,     1,
       3,     1,     3,     6,     6,     6,     8,     6,     6,     8,
      13,    13,    13,    15,    15,    15,    15,    17,     1,     1,
       4,     0,     2,     1,     3,     0,     1,     3,     1,     1,
       4,     0,     2,     1,     3,     3,     1,     1,     4,     0,
       2,     1,     3,     3,     1,     1,     4,     0,     2,     1,
       3,     3,     0,     2,     0,     2,     0,     2,     0,     2,
       6,     3,     4,     1,     3,     3,     1,     4,     1,     4,
       3,     3,     0,     2,     1,     3,     0,     2,     1,     3,
       0,     2,     1,     3,     0,     2,     1,     3,     1,     3,
       1,     3,     0,     2,     1,     1,     1,     1,     1,     4,
       0,     2,     1,     3,     1,     4,     1,     1,     0,     1,
       2,     3,     4,     1,     1,     3,     1,     3,     1,     1,
       1,     1,     1,     4,     1
};


#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)
#define YYEMPTY         (-2)
#define YYEOF           0

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                  \
do                                                              \
  if (yychar == YYEMPTY)                                        \
    {                                                           \
      yychar = (Token);                                         \
      yylval = (Value);                                         \
      YYPOPSTACK (yylen);                                       \
      yystate = *yyssp;                                         \
      goto yybackup;                                            \
    }                                                           \
  else                                                          \
    {                                                           \
      yyerror (parm, YY_("syntax error: cannot back up")); \
      YYERROR;                                                  \
    }                                                           \
while (0)

/* Error token number */
#define YYTERROR        1
#define YYERRCODE       256



/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)

/* This macro is provided for backward compatibility. */
#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
#endif


# define YY_SYMBOL_PRINT(Title, Type, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Type, Value, parm); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*----------------------------------------.
| Print this symbol's value on YYOUTPUT.  |
`----------------------------------------*/

static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, void *parm)
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  YYUSE (parm);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
  YYUSE (yytype);
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, void *parm)
{
  YYFPRINTF (yyoutput, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep, parm);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yytype_int16 *yyssp, YYSTYPE *yyvsp, int yyrule, void *parm)
{
  unsigned long int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       yystos[yyssp[yyi + 1 - yynrhs]],
                       &(yyvsp[(yyi + 1) - (yynrhs)])
                                              , parm);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule, parm); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif


#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
yystrlen (const char *yystr)
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
yystpcpy (char *yydest, const char *yysrc)
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
        switch (*++yyp)
          {
          case '\'':
          case ',':
            goto do_not_strip_quotes;

          case '\\':
            if (*++yyp != '\\')
              goto do_not_strip_quotes;
            /* Fall through.  */
          default:
            if (yyres)
              yyres[yyn] = *yyp;
            yyn++;
            break;

          case '"':
            if (yyres)
              yyres[yyn] = '\0';
            return yyn;
          }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULLPTR, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULLPTR;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULLPTR, yytname[yyx]);
                  if (! (yysize <= yysize1
                         && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                    return 2;
                  yysize = yysize1;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T yysize1 = yysize + yystrlen (yyformat);
    if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
      return 2;
    yysize = yysize1;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep, void *parm)
{
  YYUSE (yyvaluep);
  YYUSE (parm);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yytype);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}




/*----------.
| yyparse.  |
`----------*/

int
yyparse (void *parm)
{
/* The lookahead symbol.  */
int yychar;


/* The semantic value of the lookahead symbol.  */
/* Default value used for initialization, for pacifying older GCCs
   or non-GCC compilers.  */
YY_INITIAL_VALUE (static YYSTYPE yyval_default;)
YYSTYPE yylval YY_INITIAL_VALUE (= yyval_default);

    /* Number of syntax errors so far.  */
    int yynerrs;

    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       'yyss': related to states.
       'yyvs': related to semantic values.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        YYSTYPE *yyvs1 = yyvs;
        yytype_int16 *yyss1 = yyss;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * sizeof (*yyssp),
                    &yyvs1, yysize * sizeof (*yyvsp),
                    &yystacksize);

        yyss = yyss1;
        yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yytype_int16 *yyss1 = yyss;
        union yyalloc *yyptr =
          (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
        if (! yyptr)
          goto yyexhaustedlab;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
                  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      //yychar = yylex (&yylval);
     //This is taken from the code generated by Bison 2.5
     yychar = yylex (&yylval, YYLEX_PARAM);
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 7:
#line 448 "parser.yxx" /* yacc.c:1646  */
    { initfg(static_cast<ParserState*>(parm)); }
#line 1810 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 8:
#line 450 "parser.yxx" /* yacc.c:1646  */
    { initfg(static_cast<ParserState*>(parm)); }
#line 1816 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 33:
#line 508 "parser.yxx" /* yacc.c:1646  */
    {
        ParserState* pp = static_cast<ParserState*>(parm);
        yyassert(pp, !(yyvsp[-4].oSet)() || !(yyvsp[-4].oSet).some()->empty(), "Empty var int domain.");
        bool print = (yyvsp[-1].argVec)->hasAtom("output_var");
        bool introduced = (yyvsp[-1].argVec)->hasAtom("var_is_introduced");
        pp->intvarTable.put((yyvsp[-2].sValue), pp->intvars.size());
        if (print) {
          pp->output(std::string((yyvsp[-2].sValue)), new AST::IntVar(pp->intvars.size()));
        }
        if ((yyvsp[0].oArg)()) {
          AST::Node* arg = (yyvsp[0].oArg).some();
          if (arg->isInt()) {
            pp->intvars.push_back(varspec((yyvsp[-2].sValue),
              new IntVarSpec(arg->getInt(),introduced)));
          } else if (arg->isIntVar()) {
            pp->intvars.push_back(varspec((yyvsp[-2].sValue),
              new IntVarSpec(Alias(arg->getIntVar()),introduced)));
          } else {
            yyassert(pp, false, "Invalid var int initializer.");
          }
          if (!pp->hadError)
            addDomainConstraint(pp, "int_in",
                                new AST::IntVar(pp->intvars.size()-1), (yyvsp[-4].oSet));
          delete arg;
        } else {
          pp->intvars.push_back(varspec((yyvsp[-2].sValue), new IntVarSpec((yyvsp[-4].oSet),introduced)));
        }
        delete (yyvsp[-1].argVec); free((yyvsp[-2].sValue));
      }
#line 1850 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 34:
#line 538 "parser.yxx" /* yacc.c:1646  */
    {
        ParserState* pp = static_cast<ParserState*>(parm);
        bool print = (yyvsp[-1].argVec)->hasAtom("output_var");
        bool introduced = (yyvsp[-1].argVec)->hasAtom("var_is_introduced");
        pp->boolvarTable.put((yyvsp[-2].sValue), pp->boolvars.size());
        if (print) {
          pp->output(std::string((yyvsp[-2].sValue)), new AST::BoolVar(pp->boolvars.size()));
        }
        if ((yyvsp[0].oArg)()) {
          AST::Node* arg = (yyvsp[0].oArg).some();
          if (arg->isBool()) {
            pp->boolvars.push_back(varspec((yyvsp[-2].sValue),
              new BoolVarSpec(arg->getBool(),introduced)));
          } else if (arg->isBoolVar()) {
            pp->boolvars.push_back(varspec((yyvsp[-2].sValue),
              new BoolVarSpec(Alias(arg->getBoolVar()),introduced)));
          } else {
            yyassert(pp, false, "Invalid var bool initializer.");
          }
          if (!pp->hadError)
            addDomainConstraint(pp, "int_in",
                                new AST::BoolVar(pp->boolvars.size()-1), (yyvsp[-4].oSet));
          delete arg;
        } else {
          pp->boolvars.push_back(varspec((yyvsp[-2].sValue), new BoolVarSpec((yyvsp[-4].oSet),introduced)));
        }
        delete (yyvsp[-1].argVec); free((yyvsp[-2].sValue));
      }
#line 1883 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 35:
#line 567 "parser.yxx" /* yacc.c:1646  */
    { ParserState* pp = static_cast<ParserState*>(parm);
        yyassert(pp, false, "Floats not supported.");
        delete (yyvsp[-1].argVec); free((yyvsp[-2].sValue));
      }
#line 1892 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 36:
#line 572 "parser.yxx" /* yacc.c:1646  */
    {
        ParserState* pp = static_cast<ParserState*>(parm);
        bool print = (yyvsp[-1].argVec)->hasAtom("output_var");
        bool introduced = (yyvsp[-1].argVec)->hasAtom("var_is_introduced");
        pp->setvarTable.put((yyvsp[-2].sValue), pp->setvars.size());
        if (print) {
          pp->output(std::string((yyvsp[-2].sValue)), new AST::SetVar(pp->setvars.size()));
        }
        if ((yyvsp[0].oArg)()) {
          AST::Node* arg = (yyvsp[0].oArg).some();
          if (arg->isSet()) {
            pp->setvars.push_back(varspec((yyvsp[-2].sValue),
              new SetVarSpec(arg->getSet(),introduced)));
          } else if (arg->isSetVar()) {
            pp->setvars.push_back(varspec((yyvsp[-2].sValue),
              new SetVarSpec(Alias(arg->getSetVar()),introduced)));
            delete arg;
          } else {
            yyassert(pp, false, "Invalid var set initializer.");
            delete arg;
          }
          if (!pp->hadError)
            addDomainConstraint(pp, "set_subset",
                                new AST::SetVar(pp->setvars.size()-1), (yyvsp[-4].oSet));
        } else {
          pp->setvars.push_back(varspec((yyvsp[-2].sValue), new SetVarSpec((yyvsp[-4].oSet),introduced)));
        }
        delete (yyvsp[-1].argVec); free((yyvsp[-2].sValue));
      }
#line 1926 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 37:
#line 602 "parser.yxx" /* yacc.c:1646  */
    {
        ParserState* pp = static_cast<ParserState*>(parm);
        yyassert(pp, !(yyvsp[-5].oSet)() || !(yyvsp[-5].oSet).some()->empty(), "Empty int domain.");
        yyassert(pp, (yyvsp[0].arg)->isInt(), "Invalid int initializer.");
        int i = -1;
        bool isInt = (yyvsp[0].arg)->isInt(i);
        if ((yyvsp[-5].oSet)() && isInt) {
          AST::SetLit* sl = (yyvsp[-5].oSet).some();
          if (sl->interval) {
            yyassert(pp, i >= sl->min && i <= sl->max, "Empty int domain.");
          } else {
            bool found = false;
            for (unsigned int j=0; j<sl->s.size(); j++)
              if (sl->s[j] == i) {
                found = true;
                break;
              }
            yyassert(pp, found, "Empty int domain.");
          }
        }
        pp->intvals.put((yyvsp[-3].sValue), i);
        delete (yyvsp[-2].argVec); free((yyvsp[-3].sValue));
      }
#line 1954 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 38:
#line 626 "parser.yxx" /* yacc.c:1646  */
    {
        ParserState* pp = static_cast<ParserState*>(parm);
        yyassert(pp, (yyvsp[0].arg)->isBool(), "Invalid bool initializer.");
        if ((yyvsp[0].arg)->isBool()) {
          pp->boolvals.put((yyvsp[-3].sValue), (yyvsp[0].arg)->getBool());
        }
        delete (yyvsp[-2].argVec); free((yyvsp[-3].sValue));
      }
#line 1967 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 39:
#line 635 "parser.yxx" /* yacc.c:1646  */
    {
        ParserState* pp = static_cast<ParserState*>(parm);
        yyassert(pp, !(yyvsp[-5].oSet)() || !(yyvsp[-5].oSet).some()->empty(), "Empty set domain.");
        yyassert(pp, (yyvsp[0].arg)->isSet(), "Invalid set initializer.");
        AST::SetLit* set = NULL;
        if ((yyvsp[0].arg)->isSet())
          set = (yyvsp[0].arg)->getSet();
        pp->setvals.put((yyvsp[-3].sValue), *set);
        delete set;
        delete (yyvsp[-2].argVec); free((yyvsp[-3].sValue));
      }
#line 1983 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 40:
#line 648 "parser.yxx" /* yacc.c:1646  */
    {
        ParserState* pp = static_cast<ParserState*>(parm);
        yyassert(pp, (yyvsp[-10].iValue)==1, "Arrays must start at 1");
        if (!pp->hadError) {
          bool print = (yyvsp[-1].argVec)->hasCall("output_array");
          vector<int> vars((yyvsp[-8].iValue));
          yyassert(pp, !(yyvsp[-4].oSet)() || !(yyvsp[-4].oSet).some()->empty(),
                   "Empty var int domain.");
          if (!pp->hadError) {
            if ((yyvsp[0].oVarSpecVec)()) {
              vector<VarSpec*>* vsv = (yyvsp[0].oVarSpecVec).some();
              yyassert(pp, vsv->size() == static_cast<unsigned int>((yyvsp[-8].iValue)),
                       "Initializer size does not match array dimension");
              if (!pp->hadError) {
                for (int i=0; i<(yyvsp[-8].iValue); i++) {
                  IntVarSpec* ivsv = static_cast<IntVarSpec*>((*vsv)[i]);
                  if (ivsv->alias) {
                    vars[i] = ivsv->i;
                  } else {
                    vars[i] = pp->intvars.size();
                    pp->intvars.push_back(varspec((yyvsp[-2].sValue), ivsv));
                  }
                  if (!pp->hadError && (yyvsp[-4].oSet)()) {
                    Option<AST::SetLit*> opt =
                      Option<AST::SetLit*>::some(new AST::SetLit(*(yyvsp[-4].oSet).some()));
                    addDomainConstraint(pp, "int_in",
                                        new AST::IntVar(vars[i]),
                                        opt);
                  }
                }
              }
              delete vsv;
            } else {
              IntVarSpec* ispec = new IntVarSpec((yyvsp[-4].oSet),!print);
              string arrayname = "["; arrayname += (yyvsp[-2].sValue);
              for (int i=0; i<(yyvsp[-8].iValue)-1; i++) {
                vars[i] = pp->intvars.size();
                pp->intvars.push_back(varspec(arrayname, ispec));
              }
              if((yyvsp[-8].iValue)) vars[(yyvsp[-8].iValue)-1] = pp->intvars.size();
              pp->intvars.push_back(varspec((yyvsp[-2].sValue), ispec));
            }
          }
          if (print) {
            AST::Array* a = new AST::Array();
            a->a.push_back(arrayOutput((yyvsp[-1].argVec)->getCall("output_array")));
            AST::Array* output = new AST::Array();
            for (int i=0; i<(yyvsp[-8].iValue); i++)
              output->a.push_back(new AST::IntVar(vars[i]));
            a->a.push_back(output);
            a->a.push_back(new AST::String(")"));
            pp->output(std::string((yyvsp[-2].sValue)), a);
          }
          pp->intvararrays.put((yyvsp[-2].sValue), vars);
        }
        delete (yyvsp[-1].argVec); free((yyvsp[-2].sValue));
      }
#line 2045 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 41:
#line 707 "parser.yxx" /* yacc.c:1646  */
    {
        ParserState* pp = static_cast<ParserState*>(parm);
        bool print = (yyvsp[-1].argVec)->hasCall("output_array");
        yyassert(pp, (yyvsp[-10].iValue)==1, "Arrays must start at 1");
        if (!pp->hadError) {
          vector<int> vars((yyvsp[-8].iValue));
          if ((yyvsp[0].oVarSpecVec)()) {
            vector<VarSpec*>* vsv = (yyvsp[0].oVarSpecVec).some();
            yyassert(pp, vsv->size() == static_cast<unsigned int>((yyvsp[-8].iValue)),
                     "Initializer size does not match array dimension");
            if (!pp->hadError) {
              for (int i=0; i<(yyvsp[-8].iValue); i++) {
                BoolVarSpec* bvsv = static_cast<BoolVarSpec*>((*vsv)[i]);
                if (bvsv->alias)
                  vars[i] = bvsv->i;
                else {
                  vars[i] = pp->boolvars.size();
                  pp->boolvars.push_back(varspec((yyvsp[-2].sValue), (*vsv)[i]));
                }
                if (!pp->hadError && (yyvsp[-4].oSet)()) {
                  Option<AST::SetLit*> opt =
                    Option<AST::SetLit*>::some(new AST::SetLit(*(yyvsp[-4].oSet).some()));
                  addDomainConstraint(pp, "int_in",
                                      new AST::BoolVar(vars[i]),
                                      opt);
                }
              }
            }
            delete vsv;
          } else {
            for (int i=0; i<(yyvsp[-8].iValue); i++) {
              vars[i] = pp->boolvars.size();
              pp->boolvars.push_back(varspec((yyvsp[-2].sValue),
                                             new BoolVarSpec((yyvsp[-4].oSet),!print)));
            }
          }
          if (print) {
            AST::Array* a = new AST::Array();
            a->a.push_back(arrayOutput((yyvsp[-1].argVec)->getCall("output_array")));
            AST::Array* output = new AST::Array();
            for (int i=0; i<(yyvsp[-8].iValue); i++)
              output->a.push_back(new AST::BoolVar(vars[i]));
            a->a.push_back(output);
            a->a.push_back(new AST::String(")"));
            pp->output(std::string((yyvsp[-2].sValue)), a);
          }
          pp->boolvararrays.put((yyvsp[-2].sValue), vars);
        }
        delete (yyvsp[-1].argVec); free((yyvsp[-2].sValue));
      }
#line 2100 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 42:
#line 759 "parser.yxx" /* yacc.c:1646  */
    {
        ParserState* pp = static_cast<ParserState*>(parm);
        yyassert(pp, false, "Floats not supported.");
        delete (yyvsp[-1].argVec); free((yyvsp[-2].sValue));
      }
#line 2110 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 43:
#line 766 "parser.yxx" /* yacc.c:1646  */
    {
        ParserState* pp = static_cast<ParserState*>(parm);
        bool print = (yyvsp[-1].argVec)->hasCall("output_array");
        yyassert(pp, (yyvsp[-12].iValue)==1, "Arrays must start at 1");
        if (!pp->hadError) {
          vector<int> vars((yyvsp[-10].iValue));
          if ((yyvsp[0].oVarSpecVec)()) {
            vector<VarSpec*>* vsv = (yyvsp[0].oVarSpecVec).some();
            yyassert(pp, vsv->size() == static_cast<unsigned int>((yyvsp[-10].iValue)),
                     "Initializer size does not match array dimension");
            if (!pp->hadError) {
              for (int i=0; i<(yyvsp[-10].iValue); i++) {
                SetVarSpec* svsv = static_cast<SetVarSpec*>((*vsv)[i]);
                if (svsv->alias)
                  vars[i] = svsv->i;
                else {
                  vars[i] = pp->setvars.size();
                  pp->setvars.push_back(varspec((yyvsp[-2].sValue), (*vsv)[i]));
                }
                if (!pp->hadError && (yyvsp[-4].oSet)()) {
                  Option<AST::SetLit*> opt =
                    Option<AST::SetLit*>::some(new AST::SetLit(*(yyvsp[-4].oSet).some()));
                  addDomainConstraint(pp, "set_subset",
                                      new AST::SetVar(vars[i]),
                                      opt);
                }
              }
            }
            delete vsv;
          } else {
            SetVarSpec* ispec = new SetVarSpec((yyvsp[-4].oSet),!print);
            string arrayname = "["; arrayname += (yyvsp[-2].sValue);
            for (int i=0; i<(yyvsp[-10].iValue)-1; i++) {
              vars[i] = pp->setvars.size();
              pp->setvars.push_back(varspec(arrayname, ispec));
            }
            vars[(yyvsp[-10].iValue)-1] = pp->setvars.size();
            pp->setvars.push_back(varspec((yyvsp[-2].sValue), ispec));
          }
          if (print) {
            AST::Array* a = new AST::Array();
            a->a.push_back(arrayOutput((yyvsp[-1].argVec)->getCall("output_array")));
            AST::Array* output = new AST::Array();
            for (int i=0; i<(yyvsp[-10].iValue); i++)
              output->a.push_back(new AST::SetVar(vars[i]));
            a->a.push_back(output);
            a->a.push_back(new AST::String(")"));
            pp->output(std::string((yyvsp[-2].sValue)), a);
          }
          pp->setvararrays.put((yyvsp[-2].sValue), vars);
        }
        delete (yyvsp[-1].argVec); free((yyvsp[-2].sValue));
      }
#line 2168 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 44:
#line 821 "parser.yxx" /* yacc.c:1646  */
    {
        ParserState* pp = static_cast<ParserState*>(parm);
        yyassert(pp, (yyvsp[-12].iValue)==1, "Arrays must start at 1");
        yyassert(pp, (yyvsp[-1].setValue)->size() == static_cast<unsigned int>((yyvsp[-10].iValue)),
                 "Initializer size does not match array dimension");
        if (!pp->hadError)
          pp->intvalarrays.put((yyvsp[-5].sValue), *(yyvsp[-1].setValue));
        delete (yyvsp[-1].setValue);
        free((yyvsp[-5].sValue));
        delete (yyvsp[-4].argVec);
      }
#line 2184 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 45:
#line 834 "parser.yxx" /* yacc.c:1646  */
    {
        ParserState* pp = static_cast<ParserState*>(parm);
        yyassert(pp, (yyvsp[-12].iValue)==1, "Arrays must start at 1");
        yyassert(pp, (yyvsp[-1].setValue)->size() == static_cast<unsigned int>((yyvsp[-10].iValue)),
                 "Initializer size does not match array dimension");
        if (!pp->hadError)
          pp->boolvalarrays.put((yyvsp[-5].sValue), *(yyvsp[-1].setValue));
        delete (yyvsp[-1].setValue);
        free((yyvsp[-5].sValue));
        delete (yyvsp[-4].argVec);
      }
#line 2200 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 46:
#line 847 "parser.yxx" /* yacc.c:1646  */
    {
        ParserState* pp = static_cast<ParserState*>(parm);
        yyassert(pp, false, "Floats not supported.");
        delete (yyvsp[-4].argVec); free((yyvsp[-5].sValue));
      }
#line 2210 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 47:
#line 854 "parser.yxx" /* yacc.c:1646  */
    {
        ParserState* pp = static_cast<ParserState*>(parm);
        yyassert(pp, (yyvsp[-14].iValue)==1, "Arrays must start at 1");
        yyassert(pp, (yyvsp[-1].setValueList)->size() == static_cast<unsigned int>((yyvsp[-12].iValue)),
                 "Initializer size does not match array dimension");
        if (!pp->hadError)
          pp->setvalarrays.put((yyvsp[-5].sValue), *(yyvsp[-1].setValueList));
        delete (yyvsp[-1].setValueList);
        delete (yyvsp[-4].argVec); free((yyvsp[-5].sValue));
      }
#line 2225 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 48:
#line 867 "parser.yxx" /* yacc.c:1646  */
    {
        (yyval.varSpec) = new IntVarSpec((yyvsp[0].iValue),false);
      }
#line 2233 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 49:
#line 871 "parser.yxx" /* yacc.c:1646  */
    {
        int v = 0;
        ParserState* pp = static_cast<ParserState*>(parm);
        if (pp->intvarTable.get((yyvsp[0].sValue), v))
          (yyval.varSpec) = new IntVarSpec(Alias(v),false);
        else {
          pp->err << "Error: undefined identifier " << (yyvsp[0].sValue)
                  << " in line no. "
                  << yyget_lineno(pp->yyscanner) << std::endl;
          pp->hadError = true;
          (yyval.varSpec) = new IntVarSpec(0,false); // keep things consistent
        }
        free((yyvsp[0].sValue));
      }
#line 2252 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 50:
#line 886 "parser.yxx" /* yacc.c:1646  */
    {
        vector<int> v;
        ParserState* pp = static_cast<ParserState*>(parm);
        if (pp->intvararrays.get((yyvsp[-3].sValue), v)) {
          yyassert(pp,static_cast<unsigned int>((yyvsp[-1].iValue)) > 0 &&
                      static_cast<unsigned int>((yyvsp[-1].iValue)) <= v.size(),
                   "array access out of bounds");
          if (!pp->hadError)
            (yyval.varSpec) = new IntVarSpec(Alias(v[(yyvsp[-1].iValue)-1]),false);
          else
            (yyval.varSpec) = new IntVarSpec(0,false); // keep things consistent
        } else {
          pp->err << "Error: undefined array identifier " << (yyvsp[-3].sValue)
                  << " in line no. "
                  << yyget_lineno(pp->yyscanner) << std::endl;
          pp->hadError = true;
          (yyval.varSpec) = new IntVarSpec(0,false); // keep things consistent
        }
        free((yyvsp[-3].sValue));
      }
#line 2277 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 51:
#line 909 "parser.yxx" /* yacc.c:1646  */
    { (yyval.varSpecVec) = new vector<VarSpec*>(0); }
#line 2283 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 52:
#line 911 "parser.yxx" /* yacc.c:1646  */
    { (yyval.varSpecVec) = (yyvsp[-1].varSpecVec); }
#line 2289 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 53:
#line 915 "parser.yxx" /* yacc.c:1646  */
    { (yyval.varSpecVec) = new vector<VarSpec*>(1); (*(yyval.varSpecVec))[0] = (yyvsp[0].varSpec); }
#line 2295 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 54:
#line 917 "parser.yxx" /* yacc.c:1646  */
    { (yyval.varSpecVec) = (yyvsp[-2].varSpecVec); (yyval.varSpecVec)->push_back((yyvsp[0].varSpec)); }
#line 2301 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 57:
#line 922 "parser.yxx" /* yacc.c:1646  */
    { (yyval.varSpecVec) = (yyvsp[-1].varSpecVec); }
#line 2307 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 58:
#line 926 "parser.yxx" /* yacc.c:1646  */
    { (yyval.varSpec) = new FloatVarSpec((yyvsp[0].dValue),false); }
#line 2313 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 59:
#line 928 "parser.yxx" /* yacc.c:1646  */
    {
        int v = 0;
        ParserState* pp = static_cast<ParserState*>(parm);
        if (pp->floatvarTable.get((yyvsp[0].sValue), v))
          (yyval.varSpec) = new FloatVarSpec(Alias(v),false);
        else {
          pp->err << "Error: undefined identifier " << (yyvsp[0].sValue)
                  << " in line no. "
                  << yyget_lineno(pp->yyscanner) << std::endl;
          pp->hadError = true;
          (yyval.varSpec) = new FloatVarSpec(0.0,false);
        }
        free((yyvsp[0].sValue));
      }
#line 2332 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 60:
#line 943 "parser.yxx" /* yacc.c:1646  */
    {
        vector<int> v;
        ParserState* pp = static_cast<ParserState*>(parm);
        if (pp->floatvararrays.get((yyvsp[-3].sValue), v)) {
          yyassert(pp,static_cast<unsigned int>((yyvsp[-1].iValue)) > 0 &&
                      static_cast<unsigned int>((yyvsp[-1].iValue)) <= v.size(),
                   "array access out of bounds");
          if (!pp->hadError)
            (yyval.varSpec) = new FloatVarSpec(Alias(v[(yyvsp[-1].iValue)-1]),false);
          else
            (yyval.varSpec) = new FloatVarSpec(0.0,false);
        } else {
          pp->err << "Error: undefined array identifier " << (yyvsp[-3].sValue)
                  << " in line no. "
                  << yyget_lineno(pp->yyscanner) << std::endl;
          pp->hadError = true;
          (yyval.varSpec) = new FloatVarSpec(0.0,false);
        }
        free((yyvsp[-3].sValue));
      }
#line 2357 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 61:
#line 966 "parser.yxx" /* yacc.c:1646  */
    { (yyval.varSpecVec) = new vector<VarSpec*>(0); }
#line 2363 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 62:
#line 968 "parser.yxx" /* yacc.c:1646  */
    { (yyval.varSpecVec) = (yyvsp[-1].varSpecVec); }
#line 2369 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 63:
#line 972 "parser.yxx" /* yacc.c:1646  */
    { (yyval.varSpecVec) = new vector<VarSpec*>(1); (*(yyval.varSpecVec))[0] = (yyvsp[0].varSpec); }
#line 2375 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 64:
#line 974 "parser.yxx" /* yacc.c:1646  */
    { (yyval.varSpecVec) = (yyvsp[-2].varSpecVec); (yyval.varSpecVec)->push_back((yyvsp[0].varSpec)); }
#line 2381 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 65:
#line 978 "parser.yxx" /* yacc.c:1646  */
    { (yyval.varSpecVec) = (yyvsp[-1].varSpecVec); }
#line 2387 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 66:
#line 982 "parser.yxx" /* yacc.c:1646  */
    { (yyval.varSpec) = new BoolVarSpec((yyvsp[0].iValue),false); }
#line 2393 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 67:
#line 984 "parser.yxx" /* yacc.c:1646  */
    {
        int v = 0;
        ParserState* pp = static_cast<ParserState*>(parm);
        if (pp->boolvarTable.get((yyvsp[0].sValue), v))
          (yyval.varSpec) = new BoolVarSpec(Alias(v),false);
        else {
          pp->err << "Error: undefined identifier " << (yyvsp[0].sValue)
                  << " in line no. "
                  << yyget_lineno(pp->yyscanner) << std::endl;
          pp->hadError = true;
          (yyval.varSpec) = new BoolVarSpec(false,false);
        }
        free((yyvsp[0].sValue));
      }
#line 2412 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 68:
#line 999 "parser.yxx" /* yacc.c:1646  */
    {
        vector<int> v;
        ParserState* pp = static_cast<ParserState*>(parm);
        if (pp->boolvararrays.get((yyvsp[-3].sValue), v)) {
          yyassert(pp,static_cast<unsigned int>((yyvsp[-1].iValue)) > 0 &&
                      static_cast<unsigned int>((yyvsp[-1].iValue)) <= v.size(),
                   "array access out of bounds");
          if (!pp->hadError)
            (yyval.varSpec) = new BoolVarSpec(Alias(v[(yyvsp[-1].iValue)-1]),false);
          else
            (yyval.varSpec) = new BoolVarSpec(false,false);
        } else {
          pp->err << "Error: undefined array identifier " << (yyvsp[-3].sValue)
                  << " in line no. "
                  << yyget_lineno(pp->yyscanner) << std::endl;
          pp->hadError = true;
          (yyval.varSpec) = new BoolVarSpec(false,false);
        }
        free((yyvsp[-3].sValue));
      }
#line 2437 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 69:
#line 1022 "parser.yxx" /* yacc.c:1646  */
    { (yyval.varSpecVec) = new vector<VarSpec*>(0); }
#line 2443 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 70:
#line 1024 "parser.yxx" /* yacc.c:1646  */
    { (yyval.varSpecVec) = (yyvsp[-1].varSpecVec); }
#line 2449 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 71:
#line 1028 "parser.yxx" /* yacc.c:1646  */
    { (yyval.varSpecVec) = new vector<VarSpec*>(1); (*(yyval.varSpecVec))[0] = (yyvsp[0].varSpec); }
#line 2455 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 72:
#line 1030 "parser.yxx" /* yacc.c:1646  */
    { (yyval.varSpecVec) = (yyvsp[-2].varSpecVec); (yyval.varSpecVec)->push_back((yyvsp[0].varSpec)); }
#line 2461 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 73:
#line 1032 "parser.yxx" /* yacc.c:1646  */
    { (yyval.varSpecVec) = (yyvsp[-1].varSpecVec); }
#line 2467 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 74:
#line 1036 "parser.yxx" /* yacc.c:1646  */
    { (yyval.varSpec) = new SetVarSpec(Option<AST::SetLit*>::some((yyvsp[0].setLit)),false); }
#line 2473 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 75:
#line 1038 "parser.yxx" /* yacc.c:1646  */
    {
        ParserState* pp = static_cast<ParserState*>(parm);
        int v = 0;
        if (pp->setvarTable.get((yyvsp[0].sValue), v))
          (yyval.varSpec) = new SetVarSpec(Alias(v),false);
        else {
          pp->err << "Error: undefined identifier " << (yyvsp[0].sValue)
                  << " in line no. "
                  << yyget_lineno(pp->yyscanner) << std::endl;
          pp->hadError = true;
          (yyval.varSpec) = new SetVarSpec(Alias(0),false);
        }
        free((yyvsp[0].sValue));
      }
#line 2492 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 76:
#line 1053 "parser.yxx" /* yacc.c:1646  */
    {
        vector<int> v;
        ParserState* pp = static_cast<ParserState*>(parm);
        if (pp->setvararrays.get((yyvsp[-3].sValue), v)) {
          yyassert(pp,static_cast<unsigned int>((yyvsp[-1].iValue)) > 0 &&
                      static_cast<unsigned int>((yyvsp[-1].iValue)) <= v.size(),
                   "array access out of bounds");
          if (!pp->hadError)
            (yyval.varSpec) = new SetVarSpec(Alias(v[(yyvsp[-1].iValue)-1]),false);
          else
            (yyval.varSpec) = new SetVarSpec(Alias(0),false);
        } else {
          pp->err << "Error: undefined array identifier " << (yyvsp[-3].sValue)
                  << " in line no. "
                  << yyget_lineno(pp->yyscanner) << std::endl;
          pp->hadError = true;
          (yyval.varSpec) = new SetVarSpec(Alias(0),false);
        }
        free((yyvsp[-3].sValue));
      }
#line 2517 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 77:
#line 1076 "parser.yxx" /* yacc.c:1646  */
    { (yyval.varSpecVec) = new vector<VarSpec*>(0); }
#line 2523 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 78:
#line 1078 "parser.yxx" /* yacc.c:1646  */
    { (yyval.varSpecVec) = (yyvsp[-1].varSpecVec); }
#line 2529 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 79:
#line 1082 "parser.yxx" /* yacc.c:1646  */
    { (yyval.varSpecVec) = new vector<VarSpec*>(1); (*(yyval.varSpecVec))[0] = (yyvsp[0].varSpec); }
#line 2535 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 80:
#line 1084 "parser.yxx" /* yacc.c:1646  */
    { (yyval.varSpecVec) = (yyvsp[-2].varSpecVec); (yyval.varSpecVec)->push_back((yyvsp[0].varSpec)); }
#line 2541 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 81:
#line 1087 "parser.yxx" /* yacc.c:1646  */
    { (yyval.varSpecVec) = (yyvsp[-1].varSpecVec); }
#line 2547 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 82:
#line 1091 "parser.yxx" /* yacc.c:1646  */
    { (yyval.oVarSpecVec) = Option<vector<VarSpec*>* >::none(); }
#line 2553 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 83:
#line 1093 "parser.yxx" /* yacc.c:1646  */
    { (yyval.oVarSpecVec) = Option<vector<VarSpec*>* >::some((yyvsp[0].varSpecVec)); }
#line 2559 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 84:
#line 1097 "parser.yxx" /* yacc.c:1646  */
    { (yyval.oVarSpecVec) = Option<vector<VarSpec*>* >::none(); }
#line 2565 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 85:
#line 1099 "parser.yxx" /* yacc.c:1646  */
    { (yyval.oVarSpecVec) = Option<vector<VarSpec*>* >::some((yyvsp[0].varSpecVec)); }
#line 2571 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 86:
#line 1103 "parser.yxx" /* yacc.c:1646  */
    { (yyval.oVarSpecVec) = Option<vector<VarSpec*>* >::none(); }
#line 2577 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 87:
#line 1105 "parser.yxx" /* yacc.c:1646  */
    { (yyval.oVarSpecVec) = Option<vector<VarSpec*>* >::some((yyvsp[0].varSpecVec)); }
#line 2583 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 88:
#line 1109 "parser.yxx" /* yacc.c:1646  */
    { (yyval.oVarSpecVec) = Option<vector<VarSpec*>* >::none(); }
#line 2589 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 89:
#line 1111 "parser.yxx" /* yacc.c:1646  */
    { (yyval.oVarSpecVec) = Option<vector<VarSpec*>* >::some((yyvsp[0].varSpecVec)); }
#line 2595 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 90:
#line 1115 "parser.yxx" /* yacc.c:1646  */
    {
        ConExpr c((yyvsp[-4].sValue), (yyvsp[-2].argVec));
        ParserState *pp = static_cast<ParserState*>(parm);
        if (!pp->hadError) {
          try {
            pp->fg->postConstraint(c, (yyvsp[0].argVec));
          } catch (FlatZinc::Error& e) {
            yyerror(pp, e.toString().c_str());
          }
        }
        delete (yyvsp[0].argVec); free((yyvsp[-4].sValue));
      }
#line 2612 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 91:
#line 1129 "parser.yxx" /* yacc.c:1646  */
    {
        ParserState *pp = static_cast<ParserState*>(parm);
        if (!pp->hadError) {
          try {
            pp->fg->solve((yyvsp[-1].argVec));
          } catch (FlatZinc::Error& e) {
            yyerror(pp, e.toString().c_str());
          }
        } else {
          delete (yyvsp[-1].argVec);
        }
      }
#line 2629 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 92:
#line 1142 "parser.yxx" /* yacc.c:1646  */
    {
        ParserState *pp = static_cast<ParserState*>(parm);
        if (!pp->hadError) {
          try {
            if ((yyvsp[-1].bValue))
              pp->fg->minimize((yyvsp[0].iValue),(yyvsp[-2].argVec));
            else
              pp->fg->maximize((yyvsp[0].iValue),(yyvsp[-2].argVec));
          } catch (FlatZinc::Error& e) {
            yyerror(pp, e.toString().c_str());
          }
        } else {
          delete (yyvsp[-2].argVec);
        }
      }
#line 2649 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 93:
#line 1164 "parser.yxx" /* yacc.c:1646  */
    { (yyval.oSet) = Option<AST::SetLit* >::none(); }
#line 2655 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 94:
#line 1166 "parser.yxx" /* yacc.c:1646  */
    { (yyval.oSet) = Option<AST::SetLit* >::some(new AST::SetLit(*(yyvsp[-1].setValue))); }
#line 2661 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 95:
#line 1168 "parser.yxx" /* yacc.c:1646  */
    {
        (yyval.oSet) = Option<AST::SetLit* >::some(new AST::SetLit((yyvsp[-2].iValue), (yyvsp[0].iValue)));
      }
#line 2669 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 96:
#line 1174 "parser.yxx" /* yacc.c:1646  */
    { (yyval.oSet) = Option<AST::SetLit* >::none(); }
#line 2675 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 97:
#line 1176 "parser.yxx" /* yacc.c:1646  */
    { bool haveTrue = false;
        bool haveFalse = false;
        for (int i=(yyvsp[-2].setValue)->size(); i--;) {
          haveTrue |= ((*(yyvsp[-2].setValue))[i] == 1);
          haveFalse |= ((*(yyvsp[-2].setValue))[i] == 0);
        }
        delete (yyvsp[-2].setValue);
        (yyval.oSet) = Option<AST::SetLit* >::some(
          new AST::SetLit(!haveFalse,haveTrue));
      }
#line 2690 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 100:
#line 1197 "parser.yxx" /* yacc.c:1646  */
    { (yyval.setLit) = new AST::SetLit(*(yyvsp[-1].setValue)); }
#line 2696 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 101:
#line 1199 "parser.yxx" /* yacc.c:1646  */
    { (yyval.setLit) = new AST::SetLit((yyvsp[-2].iValue), (yyvsp[0].iValue)); }
#line 2702 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 102:
#line 1205 "parser.yxx" /* yacc.c:1646  */
    { (yyval.setValue) = new vector<int>(0); }
#line 2708 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 103:
#line 1207 "parser.yxx" /* yacc.c:1646  */
    { (yyval.setValue) = (yyvsp[-1].setValue); }
#line 2714 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 104:
#line 1211 "parser.yxx" /* yacc.c:1646  */
    { (yyval.setValue) = new vector<int>(1); (*(yyval.setValue))[0] = (yyvsp[0].iValue); }
#line 2720 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 105:
#line 1213 "parser.yxx" /* yacc.c:1646  */
    { (yyval.setValue) = (yyvsp[-2].setValue); (yyval.setValue)->push_back((yyvsp[0].iValue)); }
#line 2726 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 106:
#line 1217 "parser.yxx" /* yacc.c:1646  */
    { (yyval.setValue) = new vector<int>(0); }
#line 2732 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 107:
#line 1219 "parser.yxx" /* yacc.c:1646  */
    { (yyval.setValue) = (yyvsp[-1].setValue); }
#line 2738 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 108:
#line 1223 "parser.yxx" /* yacc.c:1646  */
    { (yyval.setValue) = new vector<int>(1); (*(yyval.setValue))[0] = (yyvsp[0].iValue); }
#line 2744 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 109:
#line 1225 "parser.yxx" /* yacc.c:1646  */
    { (yyval.setValue) = (yyvsp[-2].setValue); (yyval.setValue)->push_back((yyvsp[0].iValue)); }
#line 2750 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 110:
#line 1229 "parser.yxx" /* yacc.c:1646  */
    { (yyval.floatSetValue) = new vector<double>(0); }
#line 2756 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 111:
#line 1231 "parser.yxx" /* yacc.c:1646  */
    { (yyval.floatSetValue) = (yyvsp[-1].floatSetValue); }
#line 2762 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 112:
#line 1235 "parser.yxx" /* yacc.c:1646  */
    { (yyval.floatSetValue) = new vector<double>(1); (*(yyval.floatSetValue))[0] = (yyvsp[0].dValue); }
#line 2768 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 113:
#line 1237 "parser.yxx" /* yacc.c:1646  */
    { (yyval.floatSetValue) = (yyvsp[-2].floatSetValue); (yyval.floatSetValue)->push_back((yyvsp[0].dValue)); }
#line 2774 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 114:
#line 1241 "parser.yxx" /* yacc.c:1646  */
    { (yyval.setValueList) = new vector<AST::SetLit>(0); }
#line 2780 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 115:
#line 1243 "parser.yxx" /* yacc.c:1646  */
    { (yyval.setValueList) = (yyvsp[-1].setValueList); }
#line 2786 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 116:
#line 1247 "parser.yxx" /* yacc.c:1646  */
    { (yyval.setValueList) = new vector<AST::SetLit>(1); (*(yyval.setValueList))[0] = *(yyvsp[0].setLit); delete (yyvsp[0].setLit); }
#line 2792 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 117:
#line 1249 "parser.yxx" /* yacc.c:1646  */
    { (yyval.setValueList) = (yyvsp[-2].setValueList); (yyval.setValueList)->push_back(*(yyvsp[0].setLit)); delete (yyvsp[0].setLit); }
#line 2798 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 118:
#line 1257 "parser.yxx" /* yacc.c:1646  */
    { (yyval.argVec) = new AST::Array((yyvsp[0].arg)); }
#line 2804 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 119:
#line 1259 "parser.yxx" /* yacc.c:1646  */
    { (yyval.argVec) = (yyvsp[-2].argVec); (yyval.argVec)->append((yyvsp[0].arg)); }
#line 2810 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 120:
#line 1263 "parser.yxx" /* yacc.c:1646  */
    { (yyval.arg) = (yyvsp[0].arg); }
#line 2816 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 121:
#line 1265 "parser.yxx" /* yacc.c:1646  */
    { (yyval.arg) = (yyvsp[-1].argVec); }
#line 2822 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 122:
#line 1269 "parser.yxx" /* yacc.c:1646  */
    { (yyval.oArg) = Option<AST::Node*>::none(); }
#line 2828 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 123:
#line 1271 "parser.yxx" /* yacc.c:1646  */
    { (yyval.oArg) = Option<AST::Node*>::some((yyvsp[0].arg)); }
#line 2834 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 124:
#line 1275 "parser.yxx" /* yacc.c:1646  */
    { (yyval.arg) = new AST::BoolLit((yyvsp[0].iValue)); }
#line 2840 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 125:
#line 1277 "parser.yxx" /* yacc.c:1646  */
    { (yyval.arg) = new AST::IntLit((yyvsp[0].iValue)); }
#line 2846 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 126:
#line 1279 "parser.yxx" /* yacc.c:1646  */
    { (yyval.arg) = new AST::FloatLit((yyvsp[0].dValue)); }
#line 2852 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 127:
#line 1281 "parser.yxx" /* yacc.c:1646  */
    { (yyval.arg) = (yyvsp[0].setLit); }
#line 2858 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 128:
#line 1283 "parser.yxx" /* yacc.c:1646  */
    {
        vector<int> as;
        ParserState* pp = static_cast<ParserState*>(parm);
        if (pp->intvararrays.get((yyvsp[0].sValue), as)) {
          AST::Array *ia = new AST::Array(as.size());
          for (int i=as.size(); i--;)
            ia->a[i] = new AST::IntVar(as[i]);
          (yyval.arg) = ia;
        } else if (pp->boolvararrays.get((yyvsp[0].sValue), as)) {
          AST::Array *ia = new AST::Array(as.size());
          for (int i=as.size(); i--;)
            ia->a[i] = new AST::BoolVar(as[i]);
          (yyval.arg) = ia;
        } else if (pp->setvararrays.get((yyvsp[0].sValue), as)) {
          AST::Array *ia = new AST::Array(as.size());
          for (int i=as.size(); i--;)
            ia->a[i] = new AST::SetVar(as[i]);
          (yyval.arg) = ia;
        } else {
          std::vector<int> is;
          std::vector<AST::SetLit> isS;
          int ival = 0;
          bool bval = false;
          if (pp->intvalarrays.get((yyvsp[0].sValue), is)) {
            AST::Array *v = new AST::Array(is.size());
            for (int i=is.size(); i--;)
              v->a[i] = new AST::IntLit(is[i]);
            (yyval.arg) = v;
          } else if (pp->boolvalarrays.get((yyvsp[0].sValue), is)) {
            AST::Array *v = new AST::Array(is.size());
            for (int i=is.size(); i--;)
              v->a[i] = new AST::BoolLit(is[i]);
            (yyval.arg) = v;
          } else if (pp->setvalarrays.get((yyvsp[0].sValue), isS)) {
            AST::Array *v = new AST::Array(isS.size());
            for (int i=isS.size(); i--;)
              v->a[i] = new AST::SetLit(isS[i]);
            (yyval.arg) = v;
          } else if (pp->intvals.get((yyvsp[0].sValue), ival)) {
            (yyval.arg) = new AST::IntLit(ival);
          } else if (pp->boolvals.get((yyvsp[0].sValue), bval)) {
            (yyval.arg) = new AST::BoolLit(bval);
          } else {
            (yyval.arg) = getVarRefArg(pp,(yyvsp[0].sValue));
          }
        }
        free((yyvsp[0].sValue));
      }
#line 2911 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 129:
#line 1332 "parser.yxx" /* yacc.c:1646  */
    {
        ParserState* pp = static_cast<ParserState*>(parm);
        int i = -1;
        yyassert(pp, (yyvsp[-1].arg)->isInt(i), "Non-integer array index.");
        if (!pp->hadError)
          (yyval.arg) = getArrayElement(static_cast<ParserState*>(parm),(yyvsp[-3].sValue),i);
        else
          (yyval.arg) = new AST::IntLit(0); // keep things consistent
        free((yyvsp[-3].sValue));
      }
#line 2926 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 130:
#line 1345 "parser.yxx" /* yacc.c:1646  */
    { (yyval.argVec) = new AST::Array(0); }
#line 2932 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 131:
#line 1347 "parser.yxx" /* yacc.c:1646  */
    { (yyval.argVec) = (yyvsp[-1].argVec); }
#line 2938 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 132:
#line 1351 "parser.yxx" /* yacc.c:1646  */
    { (yyval.argVec) = new AST::Array((yyvsp[0].arg)); }
#line 2944 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 133:
#line 1353 "parser.yxx" /* yacc.c:1646  */
    { (yyval.argVec) = (yyvsp[-2].argVec); (yyval.argVec)->append((yyvsp[0].arg)); }
#line 2950 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 134:
#line 1361 "parser.yxx" /* yacc.c:1646  */
    {
        ParserState *pp = static_cast<ParserState*>(parm);
        if (!pp->intvarTable.get((yyvsp[0].sValue), (yyval.iValue))) {
          pp->err << "Error: unknown integer variable " << (yyvsp[0].sValue)
                  << " in line no. "
                  << yyget_lineno(pp->yyscanner) << std::endl;
          pp->hadError = true;
        }
        free((yyvsp[0].sValue));
      }
#line 2965 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 135:
#line 1372 "parser.yxx" /* yacc.c:1646  */
    {
        vector<int> tmp;
        ParserState *pp = static_cast<ParserState*>(parm);
        if (!pp->intvararrays.get((yyvsp[-3].sValue), tmp)) {
          pp->err << "Error: unknown integer variable array " << (yyvsp[-3].sValue)
                  << " in line no. "
                  << yyget_lineno(pp->yyscanner) << std::endl;
          pp->hadError = true;
        }
        if ((yyvsp[-1].iValue) == 0 || static_cast<unsigned int>((yyvsp[-1].iValue)) > tmp.size()) {
          pp->err << "Error: array index out of bounds for array " << (yyvsp[-3].sValue)
                  << " in line no. "
                  << yyget_lineno(pp->yyscanner) << std::endl;
          pp->hadError = true;
        } else {
          (yyval.iValue) = tmp[(yyvsp[-1].iValue)-1];
        }
        free((yyvsp[-3].sValue));
      }
#line 2989 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 138:
#line 1402 "parser.yxx" /* yacc.c:1646  */
    { (yyval.argVec) = NULL; }
#line 2995 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 139:
#line 1404 "parser.yxx" /* yacc.c:1646  */
    { (yyval.argVec) = (yyvsp[0].argVec); }
#line 3001 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 140:
#line 1408 "parser.yxx" /* yacc.c:1646  */
    { (yyval.argVec) = new AST::Array((yyvsp[0].arg)); }
#line 3007 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 141:
#line 1410 "parser.yxx" /* yacc.c:1646  */
    { (yyval.argVec) = (yyvsp[-2].argVec); (yyval.argVec)->append((yyvsp[0].arg)); }
#line 3013 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 142:
#line 1414 "parser.yxx" /* yacc.c:1646  */
    {
        (yyval.arg) = new AST::Call((yyvsp[-3].sValue), AST::extractSingleton((yyvsp[-1].arg))); free((yyvsp[-3].sValue));
      }
#line 3021 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 143:
#line 1418 "parser.yxx" /* yacc.c:1646  */
    { (yyval.arg) = (yyvsp[0].arg); }
#line 3027 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 144:
#line 1422 "parser.yxx" /* yacc.c:1646  */
    { (yyval.arg) = new AST::Array((yyvsp[0].arg)); }
#line 3033 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 145:
#line 1424 "parser.yxx" /* yacc.c:1646  */
    { (yyval.arg) = (yyvsp[-2].arg); (yyval.arg)->append((yyvsp[0].arg)); }
#line 3039 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 146:
#line 1428 "parser.yxx" /* yacc.c:1646  */
    { (yyval.arg) = (yyvsp[0].arg); }
#line 3045 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 147:
#line 1430 "parser.yxx" /* yacc.c:1646  */
    { (yyval.arg) = (yyvsp[-1].arg); }
#line 3051 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 148:
#line 1434 "parser.yxx" /* yacc.c:1646  */
    { (yyval.arg) = new AST::BoolLit((yyvsp[0].iValue)); }
#line 3057 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 149:
#line 1436 "parser.yxx" /* yacc.c:1646  */
    { (yyval.arg) = new AST::IntLit((yyvsp[0].iValue)); }
#line 3063 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 150:
#line 1438 "parser.yxx" /* yacc.c:1646  */
    { (yyval.arg) = new AST::FloatLit((yyvsp[0].dValue)); }
#line 3069 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 151:
#line 1440 "parser.yxx" /* yacc.c:1646  */
    { (yyval.arg) = (yyvsp[0].setLit); }
#line 3075 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 152:
#line 1442 "parser.yxx" /* yacc.c:1646  */
    {
        vector<int> as;
        ParserState* pp = static_cast<ParserState*>(parm);
        if (pp->intvararrays.get((yyvsp[0].sValue), as)) {
          AST::Array *ia = new AST::Array(as.size());
          for (int i=as.size(); i--;)
            ia->a[i] = new AST::IntVar(as[i]);
          (yyval.arg) = ia;
        } else if (pp->boolvararrays.get((yyvsp[0].sValue), as)) {
          AST::Array *ia = new AST::Array(as.size());
          for (int i=as.size(); i--;)
            ia->a[i] = new AST::BoolVar(as[i]);
          (yyval.arg) = ia;
        } else if (pp->setvararrays.get((yyvsp[0].sValue), as)) {
          AST::Array *ia = new AST::Array(as.size());
          for (int i=as.size(); i--;)
            ia->a[i] = new AST::SetVar(as[i]);
          (yyval.arg) = ia;
        } else {
          std::vector<int> is;
          int ival = 0;
          bool bval = false;
          if (pp->intvalarrays.get((yyvsp[0].sValue), is)) {
            AST::Array *v = new AST::Array(is.size());
            for (int i=is.size(); i--;)
              v->a[i] = new AST::IntLit(is[i]);
            (yyval.arg) = v;
          } else if (pp->boolvalarrays.get((yyvsp[0].sValue), is)) {
            AST::Array *v = new AST::Array(is.size());
            for (int i=is.size(); i--;)
              v->a[i] = new AST::BoolLit(is[i]);
            (yyval.arg) = v;
          } else if (pp->intvals.get((yyvsp[0].sValue), ival)) {
            (yyval.arg) = new AST::IntLit(ival);
          } else if (pp->boolvals.get((yyvsp[0].sValue), bval)) {
            (yyval.arg) = new AST::BoolLit(bval);
          } else {
            (yyval.arg) = getVarRefArg(pp,(yyvsp[0].sValue),true);
          }
        }
        free((yyvsp[0].sValue));
      }
#line 3122 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 153:
#line 1485 "parser.yxx" /* yacc.c:1646  */
    {
        ParserState* pp = static_cast<ParserState*>(parm);
        int i = -1;
        yyassert(pp, (yyvsp[-1].arg)->isInt(i), "Non-integer array index.");
        if (!pp->hadError)
          (yyval.arg) = getArrayElement(static_cast<ParserState*>(parm),(yyvsp[-3].sValue),i);
        else
          (yyval.arg) = new AST::IntLit(0); // keep things consistent
        free((yyvsp[-3].sValue));
      }
#line 3137 "parser.tab.cpp" /* yacc.c:1646  */
    break;

  case 154:
#line 1496 "parser.yxx" /* yacc.c:1646  */
    {
        (yyval.arg) = new AST::String((yyvsp[0].sValue));
        free((yyvsp[0].sValue));
      }
#line 3146 "parser.tab.cpp" /* yacc.c:1646  */
    break;


#line 3150 "parser.tab.cpp" /* yacc.c:1646  */
      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (parm, YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (parm, yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval, parm);
          yychar = YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYTERROR;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;


      yydestruct ("Error: popping",
                  yystos[yystate], yyvsp, parm);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if !defined yyoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (parm, YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval, parm);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  yystos[*yyssp], yyvsp, parm);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  return yyresult;
}
