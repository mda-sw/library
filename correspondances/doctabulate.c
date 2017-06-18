// I N C L U D E S ///////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <ctype.h>

// D E F I N E S /////////////////////////////////////////////////////////////
 
#define MAX_STRING_LENGTH		64
#define MAX_FILE_NAME           256

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
} linked_list_t;

typedef struct
{
	dword_t**		matrix;
	dword_t			num_rows;
	dword_t			num_cols;
	size_t			spacing;
	linked_list_t		row_labels;
	linked_list_t		col_labels;
} table_t;

// G L O B A L S //////////////////////////////////////////////////////////////


// F U N C T I O N S //////////////////////////////////////////////////////////

void InitList(linked_list_t *list)
{
	list->head		    = NULL;
	list->tail		    = NULL;
	list->size          = 0;
} // end ListInit

///////////////////////////////////////////////////////////////////////////////

void AddToList(linked_list_t *list, const char *str)
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

boolean_t InTableRowLabels(table_t *table, char* str, dword_t doc_component)
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
			table->matrix[curr_index][doc_component]++;
			return TRUE;
		}
		curr_index++;
	}
	

	return FALSE;
} // end InTableRowLabels

///////////////////////////////////////////////////////////////////////////////

void DestroyList(linked_list_t *list)
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

boolean_t CrossTabulate(table_t *table, char* buffer, dword_t buffer_size, dword_t doc_component)
{

	char word[MAX_STRING_LENGTH];
	boolean_t word_found;
	dword_t i, num_words, word_num_chars;

	if (!buffer)
		return(FALSE);

	num_words = 0;
	word_found = FALSE;

	for (i = 0; i < buffer_size; i++ )
	{
		if ( !word_found )
		{
			if (isalpha(buffer[i]))
			{
				word_found = TRUE;
				word_num_chars = 0;
				word[word_num_chars++] = buffer[i];
			} // end 
		} // end if
		else
		{
			if (!isalpha(buffer[i]))// && (buffer[i] != '-'))
			{
				word_found = FALSE;
				if ( word_num_chars > 1 )
				{						
					word[word_num_chars++] = '\0';

					StrToLower(word);

					InTableRowLabels(table, word, doc_component);
	
					num_words++;

				} // end if
			} // end if
			else
			{
				word[word_num_chars++] = buffer[i];
			} // end else
		} // end else
	} // end for i

	return(TRUE);
} // end CrossTabulate

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

boolean_t ExtractSection(dword_t *start_offset, dword_t *end_offset, const char *buffer, dword_t buffer_size, dword_t component_no)
{
	dword_t i, num_chars, str_length;

	char begin[64], end[64], curr_str[64];

	boolean_t start_found;


	if (!buffer)
		return(FALSE);


	sprintf(begin, "Section %d", component_no);
	sprintf(end, "Section %d", (component_no + 1));

	str_length = (dword_t)strlen(begin);

	(*start_offset) = 0;

	start_found = FALSE;

	for (i=2; i<(buffer_size-str_length); i++)
	{
		memcpy(curr_str, &buffer[i], str_length);

		curr_str[str_length] = '\0';

		if ((strcmp(curr_str, begin) == 0) && (buffer[i-1]=='\n') && (buffer[i-2]=='\n'))
		{
			(*start_offset) = i;
			start_found = TRUE;
			break;
		}
	}

	if (!start_found)
		return(FALSE);

	num_chars = 0;
	(*end_offset) = 0;
	for (i=((*start_offset) + str_length); i<(buffer_size-str_length); i++)
	{
		memcpy(curr_str, &buffer[i], str_length);
		curr_str[str_length] = '\0';
		if ((strcmp(curr_str, end) == 0) && (buffer[i-1]=='\n') && (buffer[i-2]=='\n'))
		{
			(*end_offset) = i;
			break;
		}
		num_chars++;
	}
	
	if ((*end_offset) == 0)
	{
		if (num_chars == 0)
			return(FALSE);
		else
			(*end_offset) = buffer_size;
	}

	return(TRUE);
} // end ExtractSection

///////////////////////////////////////////////////////////////////////////////

dword_t CountSections(const char *buffer, dword_t buffer_size)
{
	dword_t i, count, str_length;

	char section[64], curr_str[64];



	if (!buffer)
		return(0);


	sprintf(section, "Section 1");
	str_length = (dword_t)strlen(section);


	count = 0;
	for (i=2; i<(buffer_size-str_length); i++)
	{
		memcpy(curr_str, &buffer[i], str_length);

		curr_str[str_length] = '\0';
		if ((strcmp(curr_str, section) == 0) && (buffer[i-1]=='\n') && (buffer[i-2]=='\n'))
		{
			count++;
			sprintf(section, "Section %d", count + 1);
			str_length = (dword_t)strlen(section);
		}
	}

	return(count);
} // end CountSections

///////////////////////////////////////////////////////////////////////////////

boolean_t ExtractPart(dword_t *start_offset, dword_t *end_offset, const char *buffer, dword_t buffer_size, dword_t component_no)
{
	dword_t i, num_chars, str_length, count;

	char curr_str[64];

	boolean_t start_found;


	if (!buffer)
		return(FALSE);


	str_length = (dword_t)strlen("Part");

	(*start_offset) = 0;

	start_found = FALSE;

	count = 0;

	for (i=2; i<(buffer_size-str_length); i++)
	{
		memcpy(curr_str, &buffer[i], str_length);
		curr_str[str_length] = '\0';
		if ((strcmp(curr_str, "Part") == 0) && (buffer[i-1]=='\n') && (buffer[i-2]=='\n'))
		{
			count++;
		}

		if (component_no == count)
		{
			(*start_offset) = i;
			start_found = TRUE;
			break;
		}

	}

	if (!start_found)
		return(FALSE);

	num_chars = 0;
	(*end_offset) = 0;
	for (i=((*start_offset) + str_length); i<(buffer_size-str_length); i++)
	{
		memcpy(curr_str, &buffer[i], str_length);
		curr_str[str_length] = '\0';
		if ((strcmp(curr_str, "Part") == 0) && (buffer[i-1]=='\n') && (buffer[i-2]=='\n'))
		{
			(*end_offset) = i;
			break;
		}
		num_chars++;
	}
	
	if ((*end_offset) == 0)
	{
		if (num_chars == 0)
			return(FALSE);
		else
			(*end_offset) = buffer_size;
	}

	return(TRUE);
} // end ExtractPart

///////////////////////////////////////////////////////////////////////////////

dword_t CountParts(const char *buffer, dword_t buffer_size)
{
	dword_t i, count, str_length;

	char curr_str[64];



	if (!buffer)
		return(0);


	str_length = (dword_t)strlen("Part");

	count = 0;
	for (i=2; i<(buffer_size-str_length); i++)
	{
		memcpy(curr_str, &buffer[i], str_length);

		curr_str[str_length] = '\0';
		if ((strcmp(curr_str, "Part") == 0) && (buffer[i-1]=='\n') && (buffer[i-2]=='\n'))
		{
			count++;
		}
	}

	return(count);
} // end CountParts

///////////////////////////////////////////////////////////////////////////////

boolean_t ExtractParagraph(dword_t *start_offset, dword_t *end_offset, const char *buffer, dword_t buffer_size, dword_t component_no)
{
	dword_t i;

	boolean_t start_found;

	dword_t count;

	if (!buffer)
		return(FALSE);

	(*start_offset) = 0;

	start_found = FALSE;

	count = 0;

	for (i=0; i<(buffer_size-1); i++)
	{
		if ((buffer[i] == '\n') && (buffer[i+1] == '\n'))
		{	
			count++;		
		}

		if (component_no == count)
		{
			start_found = TRUE;
			(*start_offset) = i;
			break;
		}
	}

	if (!start_found)
		return(FALSE);

	(*end_offset) = 0;
	for (i=((*start_offset) + 2); i<(buffer_size-1); i++)
	{
		if ((buffer[i] == '\n') && (buffer[i+1] == '\n'))
		{
			(*end_offset) = i;
			break;
		}
		
	}
	
	if ((*end_offset) == 0)
	{
		return(FALSE);
	}

	return(TRUE);
} // end ExtractParagraph

///////////////////////////////////////////////////////////////////////////////

dword_t CountParagraphs(const char *buffer, dword_t buffer_size)
{
	dword_t i, count;


	if (!buffer)
		return(FALSE);


	count = 0;

	for (i=0; i<(buffer_size-1); i++)
	{
		if ((buffer[i] == '\n') && (buffer[i+1] == '\n'))
		{	
			count++;		
		}
	}

	if (count > 0)
		count--;

	return(count);
} // end CountParagraphs

///////////////////////////////////////////////////////////////////////////////

boolean_t ExtractSentence(dword_t *start_offset, dword_t *end_offset, const char *buffer, dword_t buffer_size, dword_t component_no)
{
	dword_t i;

	boolean_t start_found;

	dword_t count;

	if (!buffer)
		return(FALSE);

	(*start_offset) = 0;

	start_found = FALSE;

	count = 1;

	for (i=0; i<buffer_size; i++)
	{
		if ((buffer[i] == '.') || (buffer[i] == ';'))
		{	
			count++;		
		}

		if (component_no == count)
		{
			start_found = TRUE;
			(*start_offset) = i + 1;
			break;
		}
	}

	if (!start_found)
		return(FALSE);

	(*end_offset) = 0;
	for (i=(*start_offset); i<buffer_size; i++)
	{
		if ((buffer[i] == '.') || (buffer[i] == ';'))
		{
			(*end_offset) = i;
			break;
		}
		
	}
	
	if ((*end_offset) == 0)
	{
		return(FALSE);
	}

	return(TRUE);
} // end ExtractSentence

///////////////////////////////////////////////////////////////////////////////

boolean_t CountSentences(const char *buffer, dword_t buffer_size)
{
	dword_t i, count;


	if (!buffer)
		return(FALSE);


	count = 0;

	for (i=0; i<buffer_size; i++)
	{
		if ((buffer[i] == '.') || (buffer[i] == ';'))
		{	
			count++;		
		}
	}

	return(count);
} // end CountSentences

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
	
	char file_in_1[MAX_FILE_NAME], 
		 file_in_2[MAX_FILE_NAME], 
		 file_out_1[MAX_FILE_NAME],
		 file_out_2[MAX_FILE_NAME],
		 file_out_3[MAX_FILE_NAME],
		 file_out_4[MAX_FILE_NAME];

/*
	char file_in_1[]="arist10x.txt",
		 file_in_2[]="docwords.txt", 
		 file_out_1[]="level1.txt",
		 file_out_2[]="level2.txt",
		 file_out_3[]="level3.txt",
		 file_out_4[]="level4.txt";
*/

	linked_list_t word_list;

	node_t *p;

	table_t table;

	FILE *fp;

	char word[MAX_STRING_LENGTH], 
		 print_buffer[64], 
		 *file_buffer, 
		 *section_buffer, 
		 *part_buffer, 
		 *paragraph_buffer, 
		 *sentence_buffer;

	dword_t i, i2, i3, i4, 
		    file_size, section_size, part_size, paragraph_size, sentence_size,
			section_no, part_no, paragraph_no, sentence_no,
		    frequency,
			num_sections, num_parts, num_paragraphs, num_sentences,
		    section_start, section_end, 
			part_start, part_end,
			paragraph_start, paragraph_end,
			sentence_start, sentence_end;



    if (argv[1])
    {
        strcpy(file_in_1, argv[1]);
    }
    else
    {
        printf("parameters: document docwords [file_out]\n");
        return(0);
    }

    if (argv[2])
    {
        strcpy(file_in_2, argv[2]);
    }
    else
    {
        printf("Error: second parameter is missing\n");
        return(0);
    }

	
    if (argv[3])
	{
		sprintf(file_out_1, "%s1.txt", argv[3]);
		sprintf(file_out_2, "%s2.txt", argv[3]);
		sprintf(file_out_3, "%s3.txt", argv[3]);
		sprintf(file_out_4, "%s4.txt", argv[3]);
	}
    else
	{
		sprintf(file_out_1, "level1.txt");
		sprintf(file_out_2, "level2.txt");
		sprintf(file_out_3, "level3.txt");
		sprintf(file_out_4, "level4.txt");
	}


	// read file and store it in a buffer ///////////////////////////////////////////////

	if (!(fp = fopen(file_in_1, "rb")) )
    {
        printf("Error: cannot open input file '%s'\n", file_in_1);
		return 0;
    }

	// find file size
	fseek(fp, 0, SEEK_END);
	file_size = ftell(fp);
	fseek(fp, 0, SEEK_SET);


	// allocate memory buffer to store characters
	file_buffer = (char*)malloc(file_size);

	fread(file_buffer, sizeof(char), file_size, fp);
	fclose(fp);


	// read words and store them in a list /////////////////////////////////////////////////

	InitList(&word_list);

	if (!(fp = fopen(file_in_2, "r")) )
    {
		free(file_buffer);
        printf("Error: cannot open '%s'\n", file_in_2);
		return 0;
    }

	while(fscanf(fp, "%s %d", word, &frequency) != EOF)
	{
		AddToList(&word_list, word);
	} // end while

	fclose(fp);

	
	// create level 1 table //////////////////////////////////

	InitTable(&table);

	num_sections = CountSections(file_buffer, file_size);

	for (p = word_list.head; p != NULL; p = p->next)
	{
		AddToList(&table.row_labels, p->str);
	} // end for

	for (i=1; i<=num_sections; i++)
	{
		sprintf(print_buffer, "S%d", i);
		AddToList(&table.col_labels, print_buffer);
	}

	CreateTable(&table);

	section_no = 0;
	for (i=1; i<=num_sections; i++)
	{
		if (ExtractSection(&section_start, &section_end, file_buffer, file_size, i))
		{
			section_buffer = &file_buffer[section_start];

			section_size = section_end - section_start;

			CrossTabulate(&table, section_buffer, section_size, section_no);
			section_no++;
		} // end if
	} // end for i

	PrintToFileTable(&table, file_out_1);

	DestroyTable(&table);


	// create level 2 table //////////////////////////////////

	InitTable(&table);

	for (p = word_list.head; p != NULL; p = p->next)
	{
		AddToList(&table.row_labels, p->str);
	} // end for


	num_sections = CountSections(file_buffer, file_size);

	for (i=1; i<=num_sections; i++)
	{
		if (ExtractSection(&section_start, &section_end, file_buffer, file_size, i))
		{
			section_buffer = &file_buffer[section_start];

			section_size = section_end - section_start;

			num_parts = CountParts(section_buffer, section_size);

			for (i2=1; i2<=num_parts; i2++)
			{
				sprintf(print_buffer, "S%d%P%d", i, i2);
				AddToList(&table.col_labels, print_buffer);
			} // end for i2
		} // end if
	} // end for i


	CreateTable(&table);

	part_no = 0;
	for (i=1; i<=num_sections; i++)
	{
		if (ExtractSection(&section_start, &section_end, file_buffer, file_size, i))
		{

			section_buffer = &file_buffer[section_start];

			section_size = section_end - section_start;

			num_parts = CountParts(section_buffer, section_size);

			for (i2=1; i2<=num_parts; i2++)
			{
				if (ExtractPart(&part_start, &part_end, section_buffer, section_size, i2))
				{
					part_buffer = &section_buffer[part_start];

					part_size = part_end - part_start;

					CrossTabulate(&table, part_buffer, part_size, part_no);
					part_no++;
				} // end if
			} // end for i2
		} // end if
	} // end for i


	PrintToFileTable(&table, file_out_2);

	DestroyTable(&table);


	// create level 3 table //////////////////////////////////

	InitTable(&table);

	for (p = word_list.head; p != NULL; p = p->next)
	{
		AddToList(&table.row_labels, p->str);
	} // end for


	num_sections = CountSections(file_buffer, file_size);

	for (i=1; i<=num_sections; i++)
	{
		if (ExtractSection(&section_start, &section_end, file_buffer, file_size, i))
		{
			section_buffer = &file_buffer[section_start];

			section_size = section_end - section_start;

			num_parts = CountParts(section_buffer, section_size);

			for (i2=1; i2<=num_parts; i2++)
			{
				if (ExtractPart(&part_start, &part_end, section_buffer, section_size, i2))
				{
					part_buffer = &section_buffer[part_start];

					part_size = part_end - part_start;

					num_paragraphs = CountParagraphs(part_buffer, part_size);

					for (i3=1; i3<=num_paragraphs; i3++)
					{
						sprintf(print_buffer, "S%d%P%dA%d", i, i2, i3);
						AddToList(&table.col_labels, print_buffer);
					} // end for i3
				} // end if
			} // end for i2
		} // end if
	} // end for i


	CreateTable(&table);

	paragraph_no = 0;
	for (i=1; i<=num_sections; i++)
	{
		if (ExtractSection(&section_start, &section_end, file_buffer, file_size, i))
		{
			section_buffer = &file_buffer[section_start];

			section_size = section_end - section_start;

			num_parts = CountParts(section_buffer, section_size);

			for (i2=1; i2<=num_parts; i2++)
			{
				if (ExtractPart(&part_start, &part_end, section_buffer, section_size, i2))
				{
					part_buffer = &section_buffer[part_start];

					part_size = part_end - part_start;

					num_paragraphs = CountParagraphs(part_buffer, part_size);

					for (i3=1; i3<=num_paragraphs; i3++)
					{
						if (ExtractParagraph(&paragraph_start, &paragraph_end, part_buffer, part_size, i3))
						{
							paragraph_buffer = &part_buffer[paragraph_start];

							paragraph_size = paragraph_end - paragraph_start;

							CrossTabulate(&table, paragraph_buffer, paragraph_size, paragraph_no);
							paragraph_no++;
						}
					} // end for i3
				} // end if
			} // end for i2
		} // end if
	} // end for i

	PrintToFileTable(&table, file_out_3);

	DestroyTable(&table);


	// create level 4 table //////////////////////////////////

	InitTable(&table);

	for (p = word_list.head; p != NULL; p = p->next)
	{
		AddToList(&table.row_labels, p->str);
	} // end for


	num_sections = CountSections(file_buffer, file_size);

	for (i=1; i<=num_sections; i++)
	{
		if (ExtractSection(&section_start, &section_end, file_buffer, file_size, i))
		{
			section_buffer = &file_buffer[section_start];

			section_size = section_end - section_start;

			num_parts = CountParts(section_buffer, section_size);

			for (i2=1; i2<=num_parts; i2++)
			{
				if (ExtractPart(&part_start, &part_end, section_buffer, section_size, i2))
				{
					part_buffer = &section_buffer[part_start];

					part_size = part_end - part_start;

					num_paragraphs = CountParagraphs(part_buffer, part_size);

					for (i3=1; i3<=num_paragraphs; i3++)
					{

						if (ExtractParagraph(&paragraph_start, &paragraph_end, part_buffer, part_size, i3))
						{
							paragraph_buffer = &part_buffer[paragraph_start];

							paragraph_size = paragraph_end - paragraph_start;

							num_sentences = CountSentences(paragraph_buffer, paragraph_size);

							for (i4=1; i4<=num_sentences; i4++)
							{
								sprintf(print_buffer, "S%d%P%dA%dE%d", i, i2, i3, i4);
								AddToList(&table.col_labels, print_buffer);
							} // end for i4
						} // end if
					} // end for i3
				} // end if
			} // end for i2
		} // end if
	} // end for i


	CreateTable(&table);


	sentence_no = 0;
	for (i=1; i<=num_sections; i++)
	{
		if (ExtractSection(&section_start, &section_end, file_buffer, file_size, i))
		{
			section_buffer = &file_buffer[section_start];

			section_size = section_end - section_start;

			num_parts = CountParts(section_buffer, section_size);

			for (i2=1; i2<=num_parts; i2++)
			{
				if (ExtractPart(&part_start, &part_end, section_buffer, section_size, i2))
				{
					part_buffer = &section_buffer[part_start];

					part_size = part_end - part_start;

					num_paragraphs = CountParagraphs(part_buffer, part_size);

					for (i3=1; i3<=num_paragraphs; i3++)
					{
						if (ExtractParagraph(&paragraph_start, &paragraph_end, part_buffer, part_size, i3))
						{
							paragraph_buffer = &part_buffer[paragraph_start];

							paragraph_size = paragraph_end - paragraph_start;

							num_sentences = CountSentences(paragraph_buffer, paragraph_size);

							for (i4=1; i4<=num_sentences; i4++)
							{
								if (ExtractSentence(&sentence_start, &sentence_end, paragraph_buffer, paragraph_size, i4))
								{
									sentence_buffer = &paragraph_buffer[sentence_start];

									sentence_size = sentence_end - sentence_start;

									CrossTabulate(&table, sentence_buffer, sentence_size, sentence_no);
									sentence_no++;
								} // end if

							} // end for i4
						}
					} // end for i3
				} // end if
			} // end for i2
		} // end if
	} // end for i

	PrintToFileTable(&table, file_out_4);

	DestroyTable(&table);


	// free resources

	free(file_buffer);
	DestroyList(&word_list);

	printf("'%s', '%s', '%s', and '%s' are created\n", file_out_1, file_out_2, 
		                                               file_out_3, file_out_4);

	return(0);
}  // end main
 
