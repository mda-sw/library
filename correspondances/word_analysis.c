// I N C L U D E S ///////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <ctype.h>

// D E F I N E S /////////////////////////////////////////////////////////////
 
#define MAX_STRING_LENGTH		64
#define MAX_FILE_NAME       256

// T Y P E S /////////////////////////////////////////////////////////////////

typedef enum { FALSE = 0, TRUE = 1 } boolean_t;

typedef unsigned long dword_t;

// S T R U C T U R E S ///////////////////////////////////////////////////////

typedef struct node_s
{
	char	str[MAX_STRING_LENGTH];
	struct	node_s* next;
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
	dword_t			threshold;
} table_t;

// G L O B A L S //////////////////////////////////////////////////////////////


// F U N C T I O N S //////////////////////////////////////////////////////////

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

boolean_t InTableRowLabels(table_t *table, char* str, dword_t curr_file)
{
	node_t *p;
	dword_t curr_index;

	if (table->row_labels.head == NULL )
		return FALSE;

	curr_index = 0;

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

void StrToLower(char* str)
{
	char *p = NULL;

	for( p = str; p < str + strlen(str); p++ )
	{
		*p = tolower(*p);
	}
	
} // end StrToLower

///////////////////////////////////////////////////////////////////////////////

boolean_t GetWords(table_t *table, char *filename, dword_t curr_file)
{

	char ch, word[MAX_STRING_LENGTH];
	boolean_t word_found;
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
			if (isalpha(ch))
			{
				word_found = TRUE;
				word_num_chars = 0;
				word[word_num_chars++] = ch;
			} // end 
		} // end if
		else
		{
			if (!isalpha(ch))// && (ch != '-'))
			{
				word_found = FALSE;
				if ( word_num_chars > 1 )
				{						
					word[word_num_chars++] = '\0';

					StrToLower(word);

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

boolean_t PrintToFileTable(table_t *table, char *filename)
{
	node_t *p;
	dword_t row, col, num_files;
	FILE *fp;

	if (table->row_labels.head == NULL )
		return FALSE;

	if (!(fp = fopen(filename, "w")))
		return FALSE;


	// print colunm labels and data

	row = 0;
	for (p = table->row_labels.head; p != NULL; p = p->next)
	{
		
		num_files = 0;
		for (col = 0; col < table->num_cols; col++)
		{
			if (table->matrix[row][col])
				num_files++;
		}

		if (num_files >= table->threshold)
		{
			fprintf(fp, "%s ", p->str);

			for (col = 0; col < table->num_cols; col++)
			{
				fprintf(fp, "%d ", table->matrix[row][col]);
			}

			fprintf(fp, "%d\n", num_files);
		}

		row++;
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

//	char file_in_1[]="words.txt", 
//		 file_in_2[]="grimmstories.txt", 
//		 file_out[]="words_enhanced.txt";

	char word[MAX_STRING_LENGTH], 
		 filename[MAX_FILE_NAME];

	dword_t frequency, file_no, threshold;


    if (argv[1])
    {
        strcpy(file_in_1, argv[1]);
    }
    else
    {
        printf("\nparameters: words file_names threshold [file_out]\n");
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
    {
        threshold = atoi(argv[3]);
    }
    else
    {
        printf("Error: threshold is missing\n");
        return 0;
    }

    if (argv[4])
		strcpy(file_out, argv[4]);
    else
		strcpy(file_out, "words_enhanced.txt");


	InitTable(&table);

	table.threshold = threshold;

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
 
