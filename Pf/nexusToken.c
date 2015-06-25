#include <stdio.h>
#include "pftypes.h"
#include "nexusToken.h"

void handleComment(nexusToken *nt);
void printComment(nexusToken *nt);
void getComment(nexusToken *nt);
void skipComment(nexusToken *nt);
void getQuotedStuff(nexusToken *nt);
void getWord(nexusToken *nt);

void nextTok(nexusToken *nt)
{
	int c;

	//printf("nextTok() here.  nt->filePtr=%li, fseek=%li\n", (long int)nt->filePtr, ftell(nt->filePtr));
	//printf("start of nextTok() here.  nt->tokLen=%i, embeddedCommentLen=%i, savedCommentLen=%i\n", nt->tokLen[0], nt->embeddedCommentLen[0], nt->savedCommentLen[0]);
	//printf("c nextTok() here.  nt->max=%i\n", nt->max[0]);
	while(1) {
		c = fgetc(nt->filePtr);
		//printf("c nextTok()  got c=%c\n", c);
		if(c == EOF) {
			//printf("nextTok() got EOF\n");
			break;
		}
		if(nt->getLineEndings[0]) {
			if(c == '\n' || c == '\r') {
				nt->tokLen[0] = 1;
				nt->tok[0] = '\n';  // note only returns \n
				break;
			}
		}
		if(c == ' ' || c == '\n' || c == '\t' || c == '\r') {
			//printf("c nextTok() got space\n");
			continue;
		}
		else if(c == '[') {
			handleComment(nt);
			if(nt->tokLen[0]) {
				break;
			}
			else {
				continue;
			}
		}
		else if(c == '\'') {
			fseek(nt->filePtr, -1, 1);
			getQuotedStuff(nt);
			break;
		}
		else if(c == '(' || c == ')' || c == '[' || c == ']' || c == '{' || 
				c == '}' || c == '\\' || c == '/' || c == ',' || c == ';' || 
				c == ':' || c == '=' || c == '*' || c == '\'' || c == '"' || 
				c == '`' || c == '+' || c == '-' || c == '<' || c == '>') {
			nt->tok[0] = c;
			nt->tokLen[0] = 1;
			break;
		}
		else {
			//printf("c nextTok() about to getWord(). nt->tokLen = %i\n", nt->tokLen[0]);
			fseek(nt->filePtr, -1, 1);
			getWord(nt);
			//printf("c nextTok() after getWord(). nt->tokLen = %i\n", nt->tokLen[0]);
			break;
		}
	}
	//printf("\n");

}

void handleComment(nexusToken *nt)
{
	long int   commentStartPos;
	int       c, c2;
	int        i;
	char       p4[3];

	// We are here having read in the opening '[' character.
	commentStartPos = ftell(nt->filePtr) - 1;
	// This next would be the char following the '['
	c = fgetc(nt->filePtr);
	//printf("c handleComment() c = %c\n", c);
	if(c == EOF) {
		printf("nextToken(), handleComment().  Reached the end of the file while still in a comment.\n");
		exit(1);
	}
	if(c == ']') {
		return;
	}
	if(c == '!' && nt->writeVisibleComments[0]) {
		fseek(nt->filePtr, commentStartPos, 0);
		printComment(nt);
		return;
	}
	else if(c == '&') {
		//printf("c handleComment()  nt->getAllCommandComments[0] = %i\n", nt->getAllCommandComments[0]);
		if(nt->getAllCommandComments[0]) {
			fseek(nt->filePtr, commentStartPos, 0);
			getComment(nt);
			return;
		}
		c2 = fgetc(nt->filePtr);
		if(c2 == EOF) {
			printf("nextToken(), handleComment().  Reached the end of the file while still in a comment.\n");
			exit(1);
		}
		if(c == '&' && (c2 == 'W' || c2 == 'w') && nt->getWeightCommandComments[0]) {
			fseek(nt->filePtr, commentStartPos, 0);
			getComment(nt);
			return;
		}
		else if(c2 == '&' && nt->getP4CommandComments[0]) {
			p4[0] = 'p';
			p4[1] = '4';
			p4[2] = ' ';
			for(i = 0; i < 3; i++) {
				c2 = fgetc(nt->filePtr);
				if(c2 == EOF) {
					printf("nextToken(), handleComment().  Reached the end of the file while still in a comment.\n");
					exit(1);
				}
				if(c2 != p4[i]) {
					fseek(nt->filePtr, commentStartPos, 0);
					skipComment(nt);
					return;
				}
			}
			fseek(nt->filePtr, commentStartPos, 0);
			getComment(nt);
			return;
		}
		else {
			fseek(nt->filePtr, commentStartPos, 0);
			skipComment(nt);
			return;
		}			
	}
	else {
		fseek(nt->filePtr, commentStartPos, 0);
		skipComment(nt);
		return;
	}
	
}

void printComment(nexusToken *nt)
{
	int   level = 0;
	int   c;

	while(1) {
		c = fgetc(nt->filePtr);
		if(c == EOF) {
			printf("nextToken(), printComment().  Reached the end of the file while still in a comment.\n");
			exit(1);
		}
		printf("%c", c);
		if(c == '[') {
			level++;
		}
		else if(c == ']') {
			level--;
		}
		if(level == 0) {
			printf("\n");
			return;
		}
	}
	
}


void skipComment(nexusToken *nt)
{
	int   level = 0;
	int   c;

	while(1) {
		c = fgetc(nt->filePtr);
		if(c == EOF) {
			printf("nextToken(), skipComment().  Reached the end of the file while still in a comment.\n");
			exit(1);
		}
		if(c == '[') {
			level++;
		}
		else if(c == ']') {
			level--;
		}
		if(level == 0) {
			return;
		}
	}
}

void getComment(nexusToken *nt)
{
	int   level = 0;
	//long int startPos = ftell(nt->filePtr);
	int   c;

	if(nt->tokLen[0]) {
		while(1) {
			c = fgetc(nt->filePtr);
			if(c == EOF) {
				printf("nextToken(), getComment().  Reached the end of the file while still in a comment.\n");
				exit(1);
			}
			nt->embeddedComment[nt->embeddedCommentLen[0]] = c;
			nt->embeddedCommentLen[0]++;
			//printf("   getComment().  embeddedCommentLen[0] is %i\n", nt->embeddedCommentLen[0]);

			if(c == '[') {
				level++;
			}
			else if(c == ']') {
				level--;
			}
			if(level == 0) {
				return;
			}
		}
	}
	else {
		while(1) {
			c = fgetc(nt->filePtr);
			if(c == EOF) {
				printf("nextToken(), getComment().  Reached the end of the file while still in a comment.\n");
				exit(1);
			}
			nt->tok[nt->tokLen[0]] = c;
			nt->tokLen[0]++;

			if(c == '[') {
				level++;
			}
			else if(c == ']') {
				level--;
			}
			if(level == 0) {
				return;
			}
		}
	}
}

void getQuotedStuff(nexusToken *nt)
{
	int   c, c2;
	
	// We arrive here such that the first call to fgetc() gets the opening single quote.
	c = fgetc(nt->filePtr);
	//printf("c getQuotedStuff() c=%c tokLen=%i\n", c, nt->tokLen[0]);
	if(nt->tokLen[0]) {
		printf("nextToken(), getQuotedStuff().  tokLen should be zero.\n");
		exit(1);
	}
	nt->tok[0] = c; // Save the first char.  It will always be a single quote.
	nt->tokLen[0] = 1;

	// Check that we do not have a single quote immediately following.  2 single quotes is ok, tho.
	c = fgetc(nt->filePtr);
	if(c == EOF) {
		printf("nextToken(), getQuotedStuff().  Un-matched single quote.\n");
		printf("Reached the end of the file while still in a quote.\n");
		exit(1);
	}
	else if(c == '\'') {
		c2 = fgetc(nt->filePtr);
		//printf("c getQuotedStuff() c2=%c\n", c);
		if(c2 == EOF) {
			printf("nextToken(), getQuotedStuff().\n");
			printf("Got 2 single quotes in a row, not properly within single quotes.\n");
			printf("And then the file ended.");
			exit(1);
		}
		else if(c2 == '\'') {
			// The opening single quote has been directly followed by 2 single quotes.  Thats ok.
			fseek(nt->filePtr, -2, 1);
		}
		else {
			// The opening quote was followed by a single quote and then something else.  Thats bad.
			printf("nextToken(), getQuotedStuff().\n");
			printf("Got 2 single quotes in a row, not properly within single quotes.\n");
			exit(1);
		}
	}
	else { // The opening single quote was not followed by a single quote.  Back up one space.
		fseek(nt->filePtr, -1, 1);
	}

	// We are now at position 1 in the quote.  Read in the rest.
	while(1) {
		c = fgetc(nt->filePtr);
		//printf("c getQuotedStuff() x c=%c\n", c);
		if(c == EOF) {
			printf("nextToken(), getQuotedStuff().  Un-matched single quote.\n");
			printf("Reached the end of the file while still in a quote.\n");
			exit(1);
		}
		if(c == '\n' || c == '\r') {
			printf("nextToken(), getQuotedStuff().\n");
			printf("Got a line ending while still in a quote.\n");
			printf("This can lead to trouble, as the token can exceed the buffer, so it is not allowed.\n");
			printf("Its a bug.  Sorry.\n");
			exit(1);
		}
		else if(c == '\'') {
			nt->tok[nt->tokLen[0]] = c;
			nt->tokLen[0]++;
			c2 = fgetc(nt->filePtr);
			if(c2 == EOF) {
				break;
			}
			else if(c2 == '\'') {
				nt->tok[nt->tokLen[0]] = c2;
				nt->tokLen[0]++;
				continue;
			}
			else {
				// The char following the single quote is something
				// else.  We want to get it on the next nextTok, so
				// back up one space.
				fseek(nt->filePtr, -1, 1);
				break;
			}
		}
		else {
			nt->tok[nt->tokLen[0]] = c;
			nt->tokLen[0]++;
		}
	}
	
}

void getWord(nexusToken *nt)
{
	int   c;
	int   c2;
	int   wordIsFinished = 0;

	
	while(1) {
		c = fgetc(nt->filePtr);
		//printf("    c getWord() c=%c \n", c);
		if(c == '[') {
			//printf("c getWord() about to handleComment() \n");
			handleComment(nt);
			//printf("c getWord() returned from handleComment() \n");

			if(nt->embeddedCommentLen[0]) {
				// We got a comment, and will be returning a comment.
				// First find out whether the word is finished.			
				c2 = fgetc(nt->filePtr);
				//printf("    c2 is %c\n", c2);
				wordIsFinished = 0;
				if(c2 == EOF) {
					wordIsFinished = 1;
				}
				else if(c2 == ' ' || c2 == '\n' || c2 == '\t' || c2 == '\r' ||
						c2 == '(' || c2 == ')' || c2 == '[' || c2 == ']' || c2 == '{' || 
						c2 == '}' || c2 == '\\' || c2 == '/' || c2 == ',' || c2 == ';' || 
						c2 == ':' || c2 == '=' || c2 == '*' || c2 == '\'' || c2 == '"' || 
						c2 == '`' || c2 == '+' || c2 == '-' || c2 == '<' || c2 == '>') {
					wordIsFinished = 1;
					fseek(nt->filePtr, -1, 1); // back up 1 only if we went forward 1
				}
				else {
					fseek(nt->filePtr, -1, 1); // back up 1 only if we went forward 1
				}
				if(wordIsFinished) {
					nt->savedCommentLen[0] = nt->embeddedCommentLen[0];
					nt->embeddedCommentLen[0] = 0;
				}
				else {
					nt->savedCommentLen[0] = 0;
				}

				return;
			}
			else {
				continue;
			}
		}
		if(c == EOF) {
			break;
		}
		else if(c == ' ' || c == '\n' || c == '\t' || c == '\r' ||
				c == '(' || c == ')' || c == '[' || c == ']' || c == '{' || 
				c == '}' || c == '\\' || c == '/' || c == ',' || c == ';' || 
				c == ':' || c == '=' || c == '*' || c == '\'' || c == '"' || 
				c == '`' || c == '+' || c == '-' || c == '<' || c == '>') {
			fseek(nt->filePtr, -1, 1); // back up 1 only if we went forward 1
			break;
		}
		nt->tok[nt->tokLen[0]] = c;
		nt->tokLen[0]++;
	}
}


void nexusSkipPastNextSemiColon(nexusToken *nt)
{
	int c;

	while(1) {
		c = fgetc(nt->filePtr);
		if(c == EOF) {
			printf("nexusSkipPastNextSemiColon().\n");
			printf("Reached the end of the file without finding a semicolon.\n");
			exit(1);
		}
		else if(c == '[') {
			skipComment(nt);
		}
		else if(c == ';') {
			return;
		}
	}
}


int nexusTokenCheckLineLengths(nexusToken *nt)
{
	int c;
	int  longest = 0;
	int  thisLen = 0;

	while(1) {
		c = fgetc(nt->filePtr);
		if(c == EOF) {
			return longest;
		}
		else if(c == '\n' || c == '\r') {
			if(thisLen > longest) {
				longest = thisLen;
			}
			thisLen = 0;
		}
		else {
			thisLen++;
		}
	}
}



