#include <cstring>
#include <iostream>

struct Stack {
  char** arr = new char*[4];
  int size = 0;
  int capacity = 4;

 public:
  Stack() = default;

  ~Stack() {
    for (int i = 0; i != size; ++i) {
      delete[] arr[i];
    }
    delete[] arr;
  }
};

void IncreaseArray(char**& arr, int& capacity) {
  char** resize = new char*[capacity * 2];
  for (int i = 0; i != capacity; ++i) {
    resize[i] = arr[i];
  }
  delete[] arr;
  arr = resize;
  capacity *= 2;
}

void DecreaseArray(char**& arr, int& capacity) {
  char** resize = new char*[capacity / 2];
  for (int i = 0; i != capacity / 2; ++i) {
    resize[i] = arr[i];
  }
  delete[] arr;
  arr = resize;
  capacity /= 2;
}

void Push(Stack& stack, char* string) {
  if (stack.capacity > stack.size) {
    stack.arr[stack.size] = string;
    stack.size += 1;
    std::cout << "ok\n";
  } else {
    IncreaseArray(stack.arr, stack.capacity);
    Push(stack, string);
  }
}

void Pop(Stack& stack) {
  if (stack.size == 0) {
    std::cout << "error\n";
    return;
  }
  if (stack.size > 0.25 * stack.capacity || stack.capacity <= 4) {
    stack.size -= 1;
    std::cout << stack.arr[stack.size] << "\n";
    delete[] stack.arr[stack.size];
  } else {
    DecreaseArray(stack.arr, stack.capacity);
    Pop(stack);
  }
}

void Back(Stack& stack) {
  if (stack.size == 0) {
    std::cout << "error\n";
    return;
  }
  std::cout << stack.arr[stack.size - 1] << "\n";
}

void Size(Stack& stack) { std::cout << stack.size << "\n"; }

void Delete(Stack& stack) {
  for (int i = 0; i != stack.size; ++i) {
    delete[] stack.arr[i];
  }
  delete[] stack.arr;
}

void Clear(Stack& stack) {
  Delete(stack);
  stack.arr = new char*[4];
  stack.size = 0;
  stack.capacity = 4;
  std::cout << "ok\n";
}

char* GetCString() {
  char* temp = new char[4];
  std::memset(temp, 0, 4);
  size_t size = 0;
  size_t capacity = 3;
  int ch;
  ch = getchar();
  if (static_cast<char>(ch) == ' ') {
    ch = getchar();
  }
  while (static_cast<char>(ch) != '\n') {
    temp[size++] = static_cast<char>(ch);
    if (size == capacity) {
      char* resize = new char[(capacity + 1) * 2];
      std::copy(temp, &temp[capacity + 1], resize);
      delete[] temp;
      temp = resize;
      capacity = (capacity + 1) * 2 - 1;
    }
    ch = getchar();
  }
  temp[size] = '\0';
  return temp;
}

void AnswerQueries() {
  Stack stack;
  char string[8];
  while (true) {
    std::cin >> string;
    if (strcmp(string, "push") == 0) {
      char* temp = GetCString();
      Push(stack, temp);
    }
    if (strcmp(string, "pop") == 0) {
      Pop(stack);
    }
    if (strcmp(string, "back") == 0) {
      Back(stack);
    }
    if (strcmp(string, "size") == 0) {
      Size(stack);
    }
    if (strcmp(string, "clear") == 0) {
      Clear(stack);
    }
    if (strcmp(string, "exit") == 0) {
      std::cout << "bye\n";
      break;
    }
  }
}

int main() { AnswerQueries(); }

