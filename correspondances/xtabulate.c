// I N C L U D E S ///////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

// D E F I N E S /////////////////////////////////////////////////////////////
 
#define MAX_STRING_LENGTH		64
#define MAX_FILE_NAME           256

// Last updated at 10/05/2005 by Dimitrios Zervas

#define  FALSE 0
#define  TRUE  1

#define IS_LETTER(c) (!ispunct(c) && !isspace(c) && !isdigit(c) && (c != '-') ? 1 : 0)
  

// T Y P E S /////////////////////////////////////////////////////////////////

typedef unsigned long dword_t;

// S T R U C T U R E S ///////////////////////////////////////////////////////

typedef struct node_s
{
	char	         str[MAX_STRING_LENGTH];
	struct	node_s*  next;
} node_t;

typedef struct
{
	node_t*		head;
	node_t*		tail;
	dword_t		size;
} linked_list;

typedef struct
{
	dword_t**		matrix;
	dword_t			num_rows;
	dword_t			num_cols;
	size_t			spacing;
	linked_list		row_labels;
	linked_list		col_labels;
} table_t;

// G L O B A L S //////////////////////////////////////////////////////////////


// F U N C T I O N S //////////////////////////////////////////////////////////


void StrToLower(char* str)
{
	char *p;

	for( p = str; p < str + strlen(str); p++ )
	{
		*p = tolower(*p);
	}
	
} // end StrToLower

///////////////////////////////////////////////////////////////////////////////

void InitList(linked_list *list)
{
	list->head		    = NULL;
	list->tail		    = NULL;
	list->size          = 0;
} // end ListInit

///////////////////////////////////////////////////////////////////////////////

void AddToList(linked_list *list, const char *str)
{
	node_t *new_node;

	new_node = (node_t*)malloc(sizeof(node_t));
	
	strcpy(new_node->str, str);
	new_node->next = NULL;
	

	if (list->head == NULL)
	{
		list->head = new_node;
		list->tail = list->head;
	} // end if
	else
	{
		list->tail->next = new_node;
		list->tail = list->tail->next;
	} // end 

	list->size++;
} // end AddToList


///////////////////////////////////////////////////////////////////////////////

int InTableRowLabels(table_t *table, char* str, dword_t curr_file)
{
	node_t *p;
	dword_t curr_index;
	

	if (table->row_labels.head == NULL )
		return FALSE;

	curr_index = 0;

	StrToLower(str);
	for (p = table->row_labels.head; p != NULL; p = p->next)
	{
		if (strcmp(p->str, str) == 0)
		{
			table->matrix[curr_index][curr_file]++;
			return TRUE;
		}
		curr_index++;
	}
	

	return FALSE;
} // end InTableRowLabels

///////////////////////////////////////////////////////////////////////////////

void DestroyList(linked_list *list)
{
	node_t *p, *prev;

	p = list->head;
	while (p != NULL)
	{
		prev = p;
		p = p->next;
		if (prev)
		{
			free(prev);
			prev = NULL;
		}
	} // end while
	

	InitList(list);

} // end DestroyList



///////////////////////////////////////////////////////////////////////////////

int GetWords(table_t *table, char *filename, dword_t curr_file)
{
	char ch, word[MAX_STRING_LENGTH];
	int word_found;
	dword_t num_words, word_num_chars;
	FILE *fp;

	num_words = 0;
	word_found = FALSE;

	if (!(fp = fopen(filename, "r")))
	{
		return FALSE;
	}

	while((ch = fgetc(fp)) != EOF)
	{
		
		if ( !word_found )
		{
			if (IS_LETTER(ch))
			{
				word_found = TRUE;
				word_num_chars = 0;
				word[word_num_chars++] = ch;
			} // end 
		} // end if
		else
		{
			if (!IS_LETTER(ch))// && (ch != '-'))
			{
				word_found = FALSE;
				if ( word_num_chars > 1 )
				{						
					word[word_num_chars] = '\0';

					InTableRowLabels(table, word, curr_file);
				
					num_words++;

				} // end if
			} // end if
			else
			{
				word[word_num_chars++] = ch;
			} // end else
		} // end else
	} // end for i

	fclose(fp);

	return TRUE;

} // end GetWords

///////////////////////////////////////////////////////////////////////////////

void InitTable(table_t *table)
{
	table->num_rows = 0;
	table->num_cols = 0;

	table->matrix   = NULL;
	table->spacing  = 0;

	InitList(&table->row_labels);
	InitList(&table->col_labels);

} // end InitTable

///////////////////////////////////////////////////////////////////////////////

void CreateTable(table_t *table)
{
	dword_t row, col;

	table->num_rows = table->row_labels.size;
	table->num_cols = table->col_labels.size;

	table->matrix   = (dword_t**)malloc(sizeof(dword_t*) * table->num_rows);

	for (row=0; row<table->num_rows; row++)
	{
		table->matrix[row] = (dword_t*)malloc(sizeof(dword_t) * table->num_cols);

		for (col=0; col<table->num_cols; col++)
		{
			table->matrix[row][col] = 0;
		}
	}

} // end CreateTable

///////////////////////////////////////////////////////////////////////////////

void DestroyTable(table_t *table)
{
	dword_t row;

	if (table->matrix)
	{
		for (row=0; row<table->num_rows; row++)
		{
			if (table->matrix[row])
				free(table->matrix[row]);
		}
		free(table->matrix);
	}


	table->num_rows = 0;
	table->num_cols = 0;

	DestroyList(&table->row_labels);
	DestroyList(&table->col_labels);

	InitTable(table);

} // end CreateTable

///////////////////////////////////////////////////////////////////////////////

int PrintToFileTable(table_t *table, char *filename)
{
	node_t *p;
	dword_t row, col;
	FILE *fp;

	if (table->row_labels.head == NULL )
		return FALSE;

	if (!(fp = fopen(filename, "w")))
		return FALSE;



	// print row labels

	for (p = table->row_labels.head; p != NULL; p = p->next)
	{

		fprintf(fp, "%s ", p->str);
	}

	fprintf(fp, "\n");

	// print colunm labels and data

	col = 0;
	for (p = table->col_labels.head; p != NULL; p = p->next)
	{
		fprintf(fp, "%s ", p->str);
		
		for (row = 0; row < table->num_rows; row++)
		{
			fprintf(fp, "%d ", table->matrix[row][col]);
		
		}

		fprintf(fp, "\n");

		col++;
	}

	fclose(fp);
	
	return TRUE;
} // end PrintToFileTable

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
	table_t table;

	node_t *p;

	FILE *fp;

	char file_in_1[MAX_FILE_NAME], 
		 file_in_2[MAX_FILE_NAME], 
		 file_out[MAX_FILE_NAME];

	//char file_in_1[]="words.txt", 
	//	 file_in_2[]="grimmstories.txt", 
	//	 file_out[]="xtabulate.txt";

	char word[MAX_STRING_LENGTH], 
		 filename[MAX_FILE_NAME];

	dword_t frequency, file_no;


    if (argv[1])
    {
        strcpy(file_in_1, argv[1]);
    }
    else
    {
        printf("\nparameters: words file_names [file_out]\n");
        return 0;
    }

    if (argv[2])
    {
        strcpy(file_in_2, argv[2]);
    }
    else
    {
        printf("Error: second parameter is missing\n");
        return 0;
    }

    if (argv[3])
		strcpy(file_out, argv[3]);
    else
		strcpy(file_out, "xtabulate.txt");


	InitTable(&table);

	if (!(fp = fopen(file_in_1, "r")) )
    {
        printf("Error: cannot open %s\n", file_in_1);
		return 0;
    }

	while(fscanf(fp, "%s %d", word, &frequency) != EOF)
	{
		AddToList(&table.row_labels, word);
	} // end while

	fclose(fp);


	if (!(fp = fopen(file_in_2, "r")) )
    {
        printf("Error: cannot open %s\n", file_in_2);
		return 0;
    }
	
	while(fscanf(fp, "%s", filename) != EOF) 
	{
		AddToList(&table.col_labels, filename); 		
	}
	fclose(fp);


	CreateTable(&table);

	file_no = 0;
	for (p = table.col_labels.head; p != NULL; p = p->next)
	{
		GetWords(&table, p->str, file_no);
		file_no++;
	} // end while


	PrintToFileTable(&table, file_out);

	DestroyTable(&table);

	printf("%s is created\n", file_out);

	return 0;
}  // end main
 

