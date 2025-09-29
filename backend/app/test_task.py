from app.celery_instance import celery_app

@celery_app.task
def add(x, y):
    print(f"Received add task: {x} + {y}")
    return x + y
