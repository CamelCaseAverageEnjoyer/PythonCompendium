def foo():
    print("Сейчас __name__ =",  __name__)
    print("Работает функция foo")

if __name__ == "__main__":
    print("Сейчас __name__ =",  __name__)
    print("Работает код верхнего уровня")
    foo()
