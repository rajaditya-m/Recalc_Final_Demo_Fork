// Pretty Console.cpp : Defines the entry point for the console application.
//
#pragma once

#include <stdio.h>
#include <stdarg.h>
#include <string.h>

void PrintTextBox(char border_x, char border_y, int size_x, int size_y, const char *formatStr, ...);

void PrintTextBox(char border_x_top, char border_y_left, char border_x_bottom, char border_y_right, int size_x, int size_y, const char *formatStr, ...);

void PrintTextBox(char border_x, char border_y, const char *formatStr, ...);

void PrintTextBox(char border_x_top, char border_y_left, char border_x_bottom, char border_y_right, const char *formatStr, ...);
