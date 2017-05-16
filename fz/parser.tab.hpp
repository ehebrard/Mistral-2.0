/* A Bison parser, made by GNU Bison 2.5.  */

/* Bison interface for Yacc-like parsers in C
   
      Copyright (C) 1984, 1989-1990, 2000-2011 Free Software Foundation, Inc.
   
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


/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
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



#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{

/* Line 2068 of yacc.c  */
#line 337 "parser.yxx"
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
       


/* Line 2068 of yacc.c  */
#line 112 "parser.tab.hpp"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif




