#!/bin/bash

find . -name  '*.py' | xargs git add 
find . -name  '*.R' | xargs git add 
find . -name  '*.sh' | xargs git add 