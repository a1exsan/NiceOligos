from typing import Optional
from fastapi.responses import RedirectResponse
from nicegui import ui, app

# Временное хранилище пользователей (почта: пароль)
users = {'test@mail.ru': 'test12345'}

@ui.page('/')
def main_page() -> Optional[RedirectResponse]:
    if not app.storage.user.get('authenticated', False):
        return RedirectResponse('/login')
    with ui.column().classes('absolute-center items-center'):
        ui.label(f'Добро пожаловать, {app.storage.user["email"]}!').classes('text-2xl')
        ui.button('Выйти', on_click=logout, icon='logout').props('outline round')

def logout():
    app.storage.user.clear()
    ui.link('Войти', '/login')

@ui.page('/register')
def register_page():
    email = ui.input('Email')
    password = ui.input('Пароль', password=True)
    confirm_password = ui.input('Подтверждение пароля', password=True)

    def register():
        if not email.value or not password.value:
            ui.notify('Поля не должны быть пустыми', color='negative')
            return
        if password.value != confirm_password.value:
            ui.notify('Пароли не совпадают', color='negative')
            return
        if email.value in users:
            ui.notify('Пользователь с таким email уже зарегистрирован', color='negative')
            return
        users[email.value] = password.value
        ui.notify('Регистрация успешна!', color='positive')
        ui.link('Войти', '/login')

    ui.button('Зарегистрироваться', on_click=register)
    ui.link('Нет аккаунта? Зарегистрироваться', '/login')

@ui.page('/login')
def login_page() -> Optional[RedirectResponse]:
    if app.storage.user.get('authenticated', False):
        return RedirectResponse('/')
    email = ui.input('Email')
    password = ui.input('Пароль', password=True)

    def login():
        print(users)
        if users.get(email.value) == password.value:
            app.storage.user.update({'email': email.value, 'authenticated': True})
            ui.link('Главная', '/')
        else:
            ui.notify('Неверный email или пароль', color='negative')

    ui.button('Войти', on_click=login)
    ui.link('Нет аккаунта? Зарегистрироваться', '/register')

if __name__ in {'__main__', '__mp_main__'}:
    ui.run(storage_secret='YOUR_SECRET_KEY_HERE', port=8010)